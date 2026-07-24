#! /usr/bin/env python3
"""
Download MycoCosm (JGI fungal) genome data.

JGI retired its legacy ``get-directory`` XML download endpoint. Downloads now go
through the JGI Data Portal API (https://files.jgi.doe.gov), which this module
drives directly (no Globus):

    1. authenticate at signon.jgi.doe.gov (the ``jgi_session`` cookie value is
       the session token);
    2. list an organism's files via the ``mycocosm_file_list`` search endpoint;
    3. restore any archived (PURGED, on-tape) files via ``request_archived_files``;
    4. download immediately-available (RESTORED) files as a single zip stream via
       the ``download_files`` endpoint, authorized with
       ``Authorization: Bearer <session token>``.

The GFF3 hierarchy selects the *filtered* gene models only - JGI ``jat_label``
``genes_filtered`` (i.e. the GeneCatalog / FilteredModels ``.gff``) - never the
unfiltered ``genes_all`` models.

PLEASE respect JGI's rate limits.
"""

import os
import re
import sys
import time
import shutil
import zipfile
import logging
import argparse
import requests
import pandas as pd
from tqdm import tqdm
from urllib.parse import unquote
from mycotools.lib.kontools import format_path, outro, intro, setup_logging
from mycotools.lib.dbtools import login_check
from pathlib import Path

logger = logging.getLogger(__name__)

# JGI Data Portal API endpoints (non-Globus)
SIGNON_URL = "https://signon.jgi.doe.gov/signon/create"
SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
RESTORE_URL = "https://files.jgi.doe.gov/request_archived_files/"
DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/"


# ===========================================================================
# JGI Data Portal API - non-Globus download implementation (see module docstring)
# ===========================================================================

# Ordered, most-preferred-first jat_labels and acceptable (un-gzipped) file
# formats per download type. Most importantly, "gff3" resolves only to the
# filtered gene models (genes_filtered), never genes_all.
_TYPE_LABELS = {
    "gff3": (["genes_filtered"], {"gff", "gff3"}),
    "faa": (["proteins_filtered"], {"fasta", "fa", "aa"}),
    "transcript": (["transcripts_filtered"], {"fasta", "fa", "fna", "fsa", "nt"}),
    "est": (["ests", "est_clusters"], {"fasta", "fa", "fna", "fsa"}),
}

# Mitochondrial files that MycoCosm misfiles under a nuclear label. Mito
# assemblies are routinely tagged `assembly_unmasked` and shelved under "Genome
# Assembly (unmasked)" (e.g. Suilu4_MitoAssemblyScaffolds.fasta.gz), and mito
# annotations are tagged `genes_filtered` (e.g. Lst7536_1_MitoGenes.gff3.gz).
# Some portals list one and the same file under both the nuclear label and
# `assembly_mitochondrial`, so the label cannot discriminate and the filename
# has to. A mitochondrion is not an organismal genome and must never stand in
# for one - it is ~1/1000th the size, so the substitution silently produces a
# nonsense MTDB entry rather than an obvious failure.
#
# "mito" must open a filename token and be followed either by the rest of
# "mitochondri*" or by an assembly/annotation noun. Both guards protect nuclear
# files whose organism name merely contains the substring: it appears mid-token
# in Fomitopsis_*_AssemblyScaffolds.fasta.gz and token-initially in
# Mitosporidium_*_AssemblyScaffolds.fasta.gz.
_MITO_FILE = re.compile(
    r"(?:^|[_.\-])mito"
    r"(?:chondri\w*|(?=[_.\-]?(?:assembl|scaffold|contig|chromosom|genome|gene)))",
    re.IGNORECASE,
)


def _file_format(f):
    """Return a file's lowercase format from its metadata, falling back to the
    filename extension (ignoring a trailing .gz)."""
    fmt = ((f.get("metadata") or {}).get("file_format") or "").lower()
    if fmt:
        return fmt
    name = re.sub(r"\.gz$", "", f.get("file_name", ""), flags=re.IGNORECASE)
    ext = re.search(r"\.([^.]+)$", name)
    return ext[1].lower() if ext else ""


def _jat_label(f):
    """Return a file's lowercase JGI Analysis Task label (the canonical file
    role, e.g. assembly_masked, genes_filtered)."""
    return ((f.get("metadata") or {}).get("jat_label") or "").lower()


def _is_restored(f):
    """Whether a file is immediately downloadable (on disk) rather than PURGED
    to tape."""
    return str(f.get("file_status", "")).upper() == "RESTORED"


def _is_mito(f):
    """Whether a file is a mitochondrial assembly/annotation rather than the
    organismal one, judged by filename because JGI's labels misreport it (see
    ``_MITO_FILE``)."""
    return bool(_MITO_FILE.search(f.get("file_name", "")))


def select_file(files, ftype, masked=True):
    """Choose the single best file record for `ftype` from an organism's file
    list.

    Mitochondrial files are excluded outright, whatever their label - see
    ``_MITO_FILE``. Of the remainder, preference goes in order to:
      1. immediately-available (RESTORED) files over archived (PURGED) ones;
      2. the type's own label preference (e.g. masked assembly before unmasked
         when ``masked`` is set).

    Returns the chosen file dict, or None if the organism has no matching file.
    """
    if ftype == "fna":
        labels = (
            ["assembly_masked", "assembly_unmasked"]
            if masked
            else ["assembly_unmasked", "assembly_masked"]
        )
        formats = {"fasta", "fa", "fna", "fsa"}
        name_ok = lambda n: True
    elif ftype in _TYPE_LABELS:
        labels, formats = _TYPE_LABELS[ftype]
        if ftype == "faa":
            # proteins_filtered also tags the .tab annotation and promoter files;
            # keep only the actual proteome fasta
            name_ok = lambda n: bool(
                re.search(r"\.aa\.fa(sta)?(\.gz)?$", n, re.IGNORECASE)
            )
        else:
            name_ok = lambda n: True
    else:
        return None

    candidates = []
    for f in files:
        label = _jat_label(f)
        if label not in labels:
            continue
        if _is_mito(f):
            continue
        if _file_format(f) not in formats:
            continue
        if not name_ok(f.get("file_name", "")):
            continue
        status_rank = 0 if _is_restored(f) else 1
        candidates.append((status_rank, labels.index(label), f))
    if not candidates:
        return None
    candidates.sort(key=lambda t: (t[0], t[1]))
    return candidates[0][2]


def parse_org_name(name):
    """Split a JGI organism name (e.g. "Acaromyces ingoldii MCA 4198 v1.0") into
    (genus, species, strain), dropping a trailing version token (strain is the
    remaining words concatenated)."""
    parts = str(name).split()
    if not parts:
        return "", "", ""
    genus = parts[0]
    species = parts[1] if len(parts) > 1 else "sp."
    rest = parts[2:]
    if rest and re.fullmatch(r"[vV]?\d+(\.\d+)*", rest[-1]):
        rest = rest[:-1]
    return genus, species, "".join(rest)


def _fill_if_empty(df, i, col, value):
    """Set df.at[i, col] = value only when there is no existing non-empty value,
    so curated reference genus/species/strain are preserved while bare-accession
    input is populated from JGI metadata."""
    if not value:
        return
    cur = df.at[i, col] if col in df.columns else None
    if (
        cur is None
        or (isinstance(cur, float) and pd.isna(cur))
        or str(cur).strip() == ""
    ):
        df.at[i, col] = value


def jgi_api_login(user, pwd, max_attempts=5, spacer="\t"):
    """Authenticate against JGI's signon service and return (session, token).
    The session token (the jgi_session cookie value) authorizes the search,
    restore, and download endpoints. Exits (100) after repeated failures."""
    session = requests.Session()
    for attempt in range(1, max_attempts + 1):
        try:
            resp = session.post(
                SIGNON_URL, data={"login": user, "password": pwd}, timeout=120
            )
        except requests.RequestException as error:
            logger.warning(f"{spacer}\tJGI login error (attempt {attempt}): {error}")
            time.sleep(5)
            continue
        token = session.cookies.get("jgi_session")
        if resp.status_code == 200 and token:
            return session, unquote(token)
        logger.warning(
            f"{spacer}\tJGI login failed (attempt {attempt}, status {resp.status_code})"
        )
        time.sleep(5)
    logger.error(f"{spacer}Failed {max_attempts} JGI login attempts.")
    sys.exit(100)


def search_organism(session, portal_id, spacer="\t", max_attempts=3):
    """Return (organism_record, files) for a MycoCosm portal id (e.g. "Acain1")
    from the JGI Data Portal search endpoint, paginating to gather every file.
    Returns (None, []) if the organism is absent or the query fails."""
    org, files, page = None, [], 1
    while True:
        params = {"organism": portal_id, "api_version": "2", "x": "50", "p": str(page)}
        data = None
        for attempt in range(1, max_attempts + 1):
            try:
                resp = session.get(
                    SEARCH_URL,
                    params=params,
                    headers={"accept": "application/json"},
                    timeout=120,
                )
            except requests.RequestException as error:
                # transient per-attempt search noise during bulk assimilation;
                # DEBUG so it is hidden unless --verbose is set
                logger.debug(
                    f"{spacer}\t{portal_id} search error (attempt {attempt}): {error}"
                )
                time.sleep(2)
                continue
            if resp.status_code != 200:
                logger.debug(f"{spacer}\t{portal_id} search HTTP {resp.status_code}")
                time.sleep(2)
                continue
            try:
                data = resp.json()
            except ValueError:
                logger.debug(f"{spacer}\t{portal_id} search returned non-JSON")
                time.sleep(2)
                continue
            break
        if data is None:
            return None, []
        organisms = data.get("organisms") or []
        if not organisms:
            break
        if org is None:
            org = organisms[0]
        files.extend(organisms[0].get("files") or [])
        total = data.get("file_total") or len(files)
        if len(files) >= total or not data.get("next_page"):
            break
        page += 1
    return org, files


def _mycocosm_ids(org_id, top_hit, portal_id, file_ids):
    """Build the request_archived_files / download_files ``ids`` payload for a
    MycoCosm organism."""
    entry = {"file_ids": list(file_ids)}
    if top_hit:
        entry["top_hit"] = top_hit
    if portal_id:
        entry["mycocosm_portal_id"] = portal_id
    return {org_id: entry}


def request_restore(session, token, ids_payload, spacer="\t"):
    """Request that archived (PURGED) files be restored to disk. Returns the
    restore request's status URL, or None on failure."""
    body = {"ids": ids_payload, "send_mail": False, "api_version": "2"}
    try:
        resp = session.post(
            RESTORE_URL,
            json=body,
            headers={
                "accept": "application/json",
                "content-type": "application/json",
                "Authorization": "Bearer " + token,
            },
            timeout=120,
        )
    except requests.RequestException as error:
        logger.warning(f"{spacer}\trestore request error: {error}")
        return None
    if resp.status_code != 200:
        logger.warning(f"{spacer}\trestore request failed (HTTP {resp.status_code})")
        return None
    try:
        return resp.json().get("request_status_url")
    except ValueError:
        return None


def _fmt_elapsed(seconds):
    """Human-friendly mm:ss / h:mm elapsed string."""
    seconds = int(seconds)
    if seconds < 3600:
        return f"{seconds // 60}m{seconds % 60:02d}s"
    return f"{seconds // 3600}h{(seconds % 3600) // 60:02d}m"


# what each JGI restore status means for a tape->disk transfer, surfaced to users
_RESTORE_STATUS_MSG = {
    "new": "request queued",
    "pending": "retrieving from tape",
    "staging": "staging to disk",
    "ready": "staged to disk",
    "expired": "restore expired",
}


def poll_restore(
    session, status_url, timeout=60, interval=30, spacer="\t", label="", heartbeat=60
):
    """Poll a tape-restore request until its files are READY (returns True) or
    the timeout / expiry is reached (returns False). A `timeout` of None polls
    until the restore resolves one way or the other, however long that takes.

    Progress is logged so the user can see the tape->disk transfer advance:
    every status transition (queued -> retrieving -> staging -> ready) is
    reported, plus a heartbeat every `heartbeat` seconds while a stage lingers."""
    if not status_url:
        return False
    tag = f"{label}: " if label else ""
    waited = 0
    last_status = None
    last_heartbeat = 0
    while timeout is None or waited <= timeout:
        status = ""
        try:
            resp = session.get(
                status_url, headers={"accept": "application/json"}, timeout=60
            )
            status = (resp.json().get("status") or "").lower()
        except (requests.RequestException, ValueError):
            pass
        if status == "ready":
            logger.info(
                f"{spacer}\t{tag}tape restore complete - files staged to disk "
                f"(waited {_fmt_elapsed(waited)})"
            )
            return True
        if status == "expired":
            logger.warning(
                f"{spacer}\t{tag}tape restore expired; a new request is needed"
            )
            return False
        # surface the transfer's progress: log each stage change, then a
        # periodic heartbeat so a long-running stage does not look hung
        detail = _RESTORE_STATUS_MSG.get(status, status or "waiting")
        if status != last_status:
            logger.info(
                f"{spacer}\t{tag}tape restore: {detail} "
                f"(elapsed {_fmt_elapsed(waited)}; disk restores usually take "
                "<1 h, up to a night)"
            )
            last_status = status
            last_heartbeat = waited
        elif heartbeat and waited - last_heartbeat >= heartbeat:
            logger.info(
                f"{spacer}\t{tag}tape restore still in progress: {detail} "
                f"(elapsed {_fmt_elapsed(waited)})"
            )
            last_heartbeat = waited
        time.sleep(interval)
        waited += interval
    logger.warning(
        f"{spacer}\t{tag}tape restore did not complete within {_fmt_elapsed(timeout)}"
    )
    return False


def download_zip(session, token, ids_payload, dest_zip, spacer="\t", max_attempts=3):
    """Download the given RESTORED files as a single zip stream. Returns True on
    success (dest_zip written), False otherwise."""
    body = {"ids": ids_payload, "api_version": "2"}
    for attempt in range(1, max_attempts + 1):
        try:
            resp = session.post(
                DOWNLOAD_URL,
                json=body,
                headers={
                    "accept": "application/json",
                    "content-type": "application/json",
                    "Authorization": "Bearer " + token,
                },
                timeout=1800,
                stream=True,
            )
        except requests.RequestException as error:
            logger.warning(f"{spacer}\tdownload error (attempt {attempt}): {error}")
            time.sleep(5)
            continue
        ctype = resp.headers.get("content-type", "")
        if resp.status_code == 200 and "zip" in ctype.lower():
            with open(dest_zip, "wb") as out:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    if chunk:
                        out.write(chunk)
            resp.close()
            return True
        resp.close()
        logger.warning(
            f"{spacer}\tdownload attempt {attempt} failed (HTTP {resp.status_code}, {ctype})"
        )
        time.sleep(5)
    return False


def extract_zip(zip_path, wanted, spacer="\t"):
    """Extract files from a JGI download archive. `wanted` maps a member's
    basename -> destination path. Returns the set of basenames extracted."""
    extracted = set()
    try:
        with zipfile.ZipFile(zip_path) as archive:
            members = {Path(m).name: m for m in archive.namelist()}
            for basename, dest in wanted.items():
                member = members.get(basename)
                if member is None:
                    continue
                Path(dest).parent.mkdir(parents=True, exist_ok=True)
                with archive.open(member) as src, open(dest, "wb") as out:
                    shutil.copyfileobj(src, out)
                extracted.add(basename)
    except zipfile.BadZipFile:
        logger.error(f"{spacer}\tcorrupt JGI download archive {zip_path}")
    return extracted


def _dwnld_org(
    session, token, ids_payload, portal_id, selected, output, tmp_dir, spacer
):
    """Download one organism's `selected` {type: file record} as a single zip and
    extract it into `output/<type>/`. Returns {type: destination path} for the
    files actually obtained (empty when the download itself failed)."""
    dest_zip = os.path.join(tmp_dir, f"{portal_id}.zip")
    if not download_zip(session, token, ids_payload, dest_zip, spacer=spacer):
        logger.warning(f"{spacer}\t{portal_id}: download failed")
        return {}

    wanted, type_dest = {}, {}
    for typ, f in selected.items():
        name = f["file_name"]
        dest = os.path.join(output, typ, name)
        wanted[name] = dest
        type_dest[typ] = (name, dest)
    extracted = extract_zip(dest_zip, wanted, spacer=spacer)
    if Path(dest_zip).is_file():
        Path(dest_zip).unlink()

    obtained = {}
    for typ, (name, dest) in type_dest.items():
        if name in extracted and Path(dest).is_file():
            obtained[typ] = dest
            logger.info(f"{spacer}\t{portal_id} {typ}: {name}")
        else:
            logger.warning(f"{spacer}\t{portal_id}: {typ} missing from archive")
    return obtained


def _flush_restores(session, token, buffer, spacer="\t"):
    """Request restores for a batch of organisms in one call - the
    request_archived_files ``ids`` payload is keyed by organism, so many
    organisms ride on a single request - and stamp each with the time its
    restore was asked for, which starts its wait clock."""
    if not buffer:
        return
    ids_payload = {}
    for entry in buffer:
        ids_payload.update(entry["ids"])
    request_restore(session, token, ids_payload, spacer=spacer)
    requested_at = time.time()
    for entry in buffer:
        entry["requested_at"] = requested_at
        entry["last_request"] = requested_at
    logger.info(
        f"{spacer}\trequested tape restores for {len(buffer)} organism(s): "
        + ", ".join(e["portal_id"] for e in buffer[:5])
        + (", ..." if len(buffer) > 5 else "")
    )
    buffer.clear()


# a restore request JGI drops or lets expire would otherwise strand an
# indefinite wait forever, so a still-pending organism is re-requested this
# often (seconds)
_REREQUEST_RESTORE = 6 * 3600


def _fmt_wait(minutes):
    """Human-friendly rendering of a wait allowance given in minutes; None is an
    unbounded wait."""
    if minutes is None:
        return "as long as it takes"
    minutes = int(minutes)
    if minutes < 60:
        return f"{minutes}m"
    if minutes % 60:
        return f"{minutes // 60}h{minutes % 60:02d}m"
    return f"{minutes // 60}h"


def circle_back(
    session,
    token,
    pending,
    df,
    output,
    tmp_dir,
    dwnlds,
    ome_set,
    deferred,
    masked=True,
    restore_wait=None,
    poll_interval=60,
    spacer="\t",
):
    """Revisit organisms whose files were left on tape during the main pass,
    downloading each as JGI stages it to disk.

    Availability is re-read from ``mycocosm_file_list`` per organism - that is
    the authoritative signal. The restore request's own status URL lags behind
    the files it restored (it still reports `pending` once they are RESTORED),
    and the generic ``/search/?datasets=`` endpoint, though it accepts many
    organisms at once, silently omits some MycoCosm portals; neither can be
    trusted here.

    One sweep of the pending organisms is made per `poll_interval` seconds. An
    organism is downloaded as soon as every file it still needs is RESTORED.
    `restore_wait` is how many minutes any one organism is waited on before it
    is given up on - deferred, not failed; None (the default) waits for as long
    as JGI takes, re-requesting a restore that has gone stale."""
    wait_seconds = None if restore_wait is None else max(0, restore_wait * 60)
    logger.info(
        f"{spacer}{len(pending)} genome(s) awaiting tape restore; checking every "
        f"{poll_interval}s, waiting {_fmt_wait(restore_wait)} per genome"
    )
    started = time.time()
    while pending:
        sweep_start = time.time()
        for portal_id in list(pending):
            entry = pending[portal_id]
            requested_at = entry.get("requested_at") or sweep_start
            waited = time.time() - requested_at

            # a dropped or expired restore request would strand an unbounded
            # wait, so reissue one that has been outstanding too long
            if time.time() - (entry.get("last_request") or requested_at) >= (
                _REREQUEST_RESTORE
            ):
                logger.info(
                    f"{spacer}\t{portal_id}: still on tape after "
                    f"{_fmt_elapsed(waited)}; re-requesting its restore"
                )
                request_restore(session, token, entry["ids"], spacer=spacer)
                entry["last_request"] = time.time()

            org, files = search_organism(session, portal_id, spacer=spacer)
            still_needed = {}
            if org:
                for typ in entry["typs"]:
                    chosen = select_file(files, typ, masked=masked)
                    if chosen is not None:
                        still_needed[typ] = chosen
                if not still_needed:
                    # JGI publishes none of these files any more - a permanent
                    # condition, so do not spend the wait allowance on it
                    logger.warning(
                        f"{spacer}\t{portal_id}: JGI no longer lists the requested "
                        "file(s)"
                    )
                    _finish_org(
                        df,
                        entry["i"],
                        portal_id,
                        dwnlds,
                        entry["dwnlded"],
                        ome_set,
                        spacer,
                    )
                    del pending[portal_id]
                    continue

            if still_needed and all(_is_restored(f) for f in still_needed.values()):
                logger.info(
                    f"{spacer}\t{portal_id}: staged to disk after "
                    f"{_fmt_elapsed(waited)}; downloading"
                )
                ids_payload = _mycocosm_ids(
                    org.get("id"),
                    (org.get("top_hit") or {}).get("_id"),
                    org.get("mycocosm_portal_id") or portal_id,
                    [f["_id"] for f in still_needed.values()],
                )
                entry["dwnlded"].update(
                    _dwnld_org(
                        session,
                        token,
                        ids_payload,
                        portal_id,
                        still_needed,
                        output,
                        tmp_dir,
                        spacer,
                    )
                )
                _finish_org(
                    df, entry["i"], portal_id, dwnlds, entry["dwnlded"], ome_set, spacer
                )
                del pending[portal_id]
            elif wait_seconds is not None and waited >= wait_seconds:
                logger.warning(
                    f"{spacer}\t{portal_id}: still on tape after "
                    f"{_fmt_elapsed(waited)}; deferring to a later run"
                )
                ome_set.add(portal_id)
                deferred.add(portal_id)
                del pending[portal_id]

            # spread a sweep's status checks over the interval so a large
            # backlog never bursts requests at JGI
            if pending:
                time.sleep(min(poll_interval / len(pending), 1))

        if pending:
            # an unbounded wait can be a long one; say what it is waiting on
            logger.info(
                f"{spacer}\t{len(pending)} genome(s) still on tape after "
                f"{_fmt_elapsed(time.time() - started)}"
            )
            remaining = poll_interval - (time.time() - sweep_start)
            if remaining > 0:
                time.sleep(remaining)


def _finish_org(df, i, portal_id, dwnlds, dwnlded, ome_set, spacer="\t"):
    """Record an organism's retrieved files as ``<type>_path`` columns and flag
    it as failed when an essential type (assembly or gff3) was not obtained -
    whether JGI had no such file, or the download/extraction did not yield it."""
    for typ, path in dwnlded.items():
        df.at[i, typ + "_path"] = path
    essential = [t for t in ("fna", "gff3") if t in dwnlds and t not in dwnlded]
    if essential:
        logger.warning(
            f"{spacer}\t{portal_id}: no {'/'.join(essential)} retrieved; excluding"
        )
        ome_set.add(portal_id)


def main(
    df,
    output,
    user,
    pwd,
    assembly=True,
    proteome=False,
    gff3=True,
    transcript=False,
    est=False,
    masked=True,
    spacer="\t",
    restore_wait=None,
    poll_interval=60,
    request_delay=3,
    ome_col=None,
    deferred=None,
    defer_tape=False,
    restore_chunk=50,
):
    """Download MycoCosm data for the JGI portal ids in `df` via the JGI Data
    Portal API (non-Globus), preserving the legacy contract: `df` gains
    ``<type>_path`` columns (e.g. fna_path, gff3_path) plus genus/species/strain,
    and the function returns (df, failed_portal_ids).

    `ome_col` names the portal id column; it defaults to ``assembly_acc``, or the
    lone column of a single-column input (MycoCosm tables label it ``portal``).

    Most MycoCosm files are archived on tape (file_status PURGED) and must be
    staged to disk before they can be downloaded. `defer_tape` chooses how that
    wait is spent:

      - False (default): each organism's restore is awaited in place, up to
        `restore_wait` minutes, before moving to the next portal id.
      - True: an organism needing tape I/O is skipped immediately - its restore
        is requested (batched `restore_chunk` organisms to a call) and it is
        revisited only after every portal id has been visited, then polled every
        `poll_interval` seconds until staged. Far faster over a large table,
        where nearly every genome needs a restore.

    Either way `restore_wait` is the maximum a single genome is waited on, in
    minutes, after which it is deferred rather than failed. None (the default)
    waits for as long as JGI takes. Pass a set as `deferred` to receive the
    portal ids given up on: unlike genuine failures (portal absent from
    MycoCosm, no such file type, corrupt download) they are pending JGI tape
    I/O, so callers should retry them rather than blacklist them."""
    if deferred is None:
        deferred = set()
    if ome_col is not None:
        if ome_col not in df.columns:
            logger.error(f"Invalid input. No {ome_col} column.")
            return df, set()
    elif "assembly_acc" in df.columns:
        ome_col = "assembly_acc"
    elif len(df.columns) == 1:
        ome_col = list(df.columns)[0]
    else:
        logger.error("Invalid input. No assembly_acc column and more than one column.")
        return df, set()

    dwnlds = []
    if assembly:
        dwnlds.append("fna")
    if proteome:
        dwnlds.append("faa")
    if gff3:
        dwnlds.append("gff3")
    if transcript:
        dwnlds.append("transcript")
    if est:
        dwnlds.append("est")

    output = str(output).rstrip("/")
    for typ in dwnlds:
        Path(os.path.join(output, typ)).mkdir(parents=True, exist_ok=True)
    tmp_dir = os.path.join(output, "jgi_zip")
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    logger.info(spacer + "Logging into JGI")
    session, token = jgi_api_login(user, pwd, spacer=spacer)

    logger.info(
        f"{spacer}Downloading {len(df)} JGI organism(s) via the JGI Data Portal API"
    )
    ome_set = set()
    pending, restore_buffer = {}, []
    for i, row in tqdm(df.iterrows(), total=len(df)):
        portal_id = row[ome_col]

        org, files = search_organism(session, portal_id, spacer=spacer)
        if not org:
            # expected during bulk assimilation; DEBUG so it is hidden unless
            # --verbose is set
            logger.debug(f"{spacer}\t{portal_id} not found in JGI MycoCosm")
            ome_set.add(portal_id)
            continue

        org_id = org.get("id")
        top_hit = (org.get("top_hit") or {}).get("_id")
        portal = org.get("mycocosm_portal_id") or portal_id

        genus, species, strain = parse_org_name(org.get("name") or "")
        _fill_if_empty(df, i, "genus", genus)
        _fill_if_empty(df, i, "species", species)
        _fill_if_empty(df, i, "strain", strain)

        selected, absent = {}, []
        for typ in dwnlds:
            chosen = select_file(files, typ, masked=masked)
            if chosen is None:
                absent.append(typ)
            else:
                selected[typ] = chosen
        if absent:
            logger.warning(f"{spacer}\t{portal_id}: JGI has no {'/'.join(absent)} file")

        # resume a previous run: any file already retrieved is kept as-is.
        # Downloads arrive gzipped, but curation decompresses in place, so an
        # unzipped copy counts too
        dwnlded = {}
        for typ, f in tuple(selected.items()):
            dest = os.path.join(output, typ, f["file_name"])
            for path in (dest, re.sub(r"\.gz$", "", dest)):
                if Path(path).is_file() and Path(path).stat().st_size > 0:
                    dwnlded[typ] = path
                    del selected[typ]
                    logger.info(
                        f"{spacer}\t{portal_id} {typ}: {Path(path).name} (preexisting)"
                    )
                    break

        if not selected:
            # nothing left to retrieve - either all preexisting, or JGI has no
            # files for any requested type
            _finish_org(df, i, portal_id, dwnlds, dwnlded, ome_set, spacer)
            continue

        # JGI keeps most files in tape archive (file_status PURGED); those must
        # be transferred to disk (RESTORED) before they can be downloaded. Report
        # which files are on tape vs already on disk, then request the restore.
        on_disk = [f for f in selected.values() if _is_restored(f)]
        on_tape = [f for f in selected.values() if not _is_restored(f)]
        if on_disk:
            logger.info(
                f"{spacer}\t{portal_id}: {len(on_disk)} file(s) already on disk: "
                + ", ".join(f.get("file_name", f["_id"]) for f in on_disk)
            )
        if on_tape:
            logger.debug(
                f"{spacer}\t{portal_id}: {len(on_tape)} file(s) are archived on TAPE and "
                "must be restored to disk before download - "
                + ", ".join(f.get("file_name", f["_id"]) for f in on_tape)
            )
            tape_ids = _mycocosm_ids(
                org_id, top_hit, portal, [f["_id"] for f in on_tape]
            )
            if defer_tape:
                # do not block the pass on JGI tape I/O: queue the restore and
                # come back to this organism once every portal has been visited
                pending[portal_id] = {
                    "i": i,
                    "portal_id": portal_id,
                    "typs": list(selected),
                    "dwnlded": dwnlded,
                    "ids": tape_ids,
                    "requested_at": None,
                }
                restore_buffer.append(pending[portal_id])
                if len(restore_buffer) >= restore_chunk:
                    _flush_restores(session, token, restore_buffer, spacer=spacer)
                continue

            status_url = request_restore(session, token, tape_ids, spacer=spacer)
            if not poll_restore(
                session,
                status_url,
                timeout=None if restore_wait is None else restore_wait * 60,
                interval=poll_interval,
                spacer=spacer,
                label=portal_id,
            ):
                logger.warning(
                    f"{spacer}\t{portal_id}: tape restore still pending; deferring this "
                    "genome (rerun later to resume once JGI has staged it to disk)"
                )
                ome_set.add(portal_id)
                deferred.add(portal_id)
                continue

        ids_payload = _mycocosm_ids(
            org_id, top_hit, portal, [f["_id"] for f in selected.values()]
        )
        dwnlded.update(
            _dwnld_org(
                session,
                token,
                ids_payload,
                portal_id,
                selected,
                output,
                tmp_dir,
                spacer,
            )
        )
        _finish_org(df, i, portal_id, dwnlds, dwnlded, ome_set, spacer)
        if request_delay:
            time.sleep(request_delay)

    # circle back to the organisms skipped for tape restores above
    _flush_restores(session, token, restore_buffer, spacer=spacer)
    if pending:
        circle_back(
            session,
            token,
            pending,
            df,
            output,
            tmp_dir,
            dwnlds,
            ome_set,
            deferred,
            masked=masked,
            restore_wait=restore_wait,
            poll_interval=poll_interval,
            spacer=spacer,
        )

    # tidy the scratch zip dir if empty
    try:
        Path(tmp_dir).rmdir()
    except OSError:
        pass

    for col in ("gff3", "faa", "fna"):
        if col in df.columns:
            del df[col]

    return df, ome_set


def cli():

    parser = argparse.ArgumentParser(
        description="Imports table/database with a JGI `assembly_acc` column (MycoCosm "
        + "portal ids) and downloads assembly, proteome, and/or gff3 via the JGI Data "
        + "Portal API. Supports rerunning/continuing previous runs in the same directory."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Genome code or table with `assembly_acc` column of JGI ome codes",
    )
    parser.add_argument(
        "-a",
        "--assembly",
        default=False,
        action="store_true",
        help="Download assembly fastas",
    )
    parser.add_argument(
        "-p",
        "--proteome",
        default=False,
        action="store_true",
        help="Download proteome fastas",
    )
    parser.add_argument(
        "-g", "--gff", default=False, action="store_true", help="Download gff3s"
    )
    parser.add_argument(
        "-t",
        "--transcript",
        default=False,
        action="store_true",
        help="Download transcripts fastas",
    )
    parser.add_argument(
        "-e", "--est", default=False, action="store_true", help="Download EST fastas"
    )
    parser.add_argument(
        "--nonmasked",
        default=False,
        action="store_true",
        help="[-a] Download nonmasked assemblies",
    )
    parser.add_argument(
        "-s",
        "--skip-tape",
        default=False,
        action="store_true",
        help="Skip genomes with tape-archived files, request their restore, and "
        + "circle back to download them once JGI stages them to disk",
    )
    parser.add_argument(
        "-w",
        "--tape-wait",
        type=int,
        default=None,
        help="Maximum minutes to wait for a single genome's tape restore before "
        + "deferring it to a later run. DEFAULT: wait indefinitely",
    )
    parser.add_argument("-o", "--output", default=str(Path.cwd()), help="Output dir")
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Report per-genome search/download diagnostics (DEBUG logging)",
    )
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    if args.nonmasked:
        args.assembly = True

    if (
        not args.assembly
        and not args.proteome
        and not args.transcript
        and not args.est
        and not args.gff
    ):
        logger.error("You must choose at least one download option.")

    ncbi_api, user, pwd = login_check(ncbi=False)

    args_dict = {
        "JGI Table": args.input,
        "Assemblies": args.assembly,
        "RepeatMasked": not args.nonmasked,
        "Proteomes": args.proteome,
        ".gff3's": args.gff,
        "Transcripts": args.transcript,
        "EST": args.est,
        "Skip tape files": args.skip_tape,
        "Max tape wait": (
            "indefinite" if args.tape_wait is None else f"{args.tape_wait} minute(s)"
        ),
    }

    start_time = intro("Download JGI files", args_dict)
    logger.warning(
        "This script does NOT account for use-restricted data. "
        + "It is user responsibility to determine use restriction status "
        + "in accord with the MycoCosm terms and conditions: "
        + "https://jgi.doe.gov/user-programs/pmo-overview/policies/legacy-data-policies/"
    )
    logger.info("")

    if Path(args.input).is_file():
        with open(args.input, "r") as raw:
            for line in raw:
                if "assembly_acc" in line.rstrip().split("\t"):
                    df = pd.read_csv(args.input, sep="\t", index_col=None)
                else:
                    df = pd.read_csv(args.input, sep="\t", header=None)
                break
    else:
        in_data = args.input.replace('"', "").replace("'", "").replace(",", " ").split()
        df = pd.DataFrame({"assembly_acc": in_data})

    output = format_path(args.output)

    jgi_df, ome_set = main(
        df,
        output,
        user,
        pwd,
        args.assembly,
        args.proteome,
        args.gff,
        args.transcript,
        args.est,
        not args.nonmasked,
        spacer="",
        restore_wait=args.tape_wait,
        defer_tape=args.skip_tape,
    )
    jgi_df = jgi_df.rename(columns={"assembly_acc": "#assembly_acc"})
    jgi_df["source"] = "jgi"
    jgi_df["restriction"] = "no"
    jgi_df.to_csv(str(Path(args.input)) + ".predb.tsv", sep="\t", index=False)

    outro(start_time)


if __name__ == "__main__":
    cli()

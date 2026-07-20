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

The file-selection hierarchy mirrors ``parse_xml`` (retained below as the
canonical reference and for backwards-compatible imports). In particular the
GFF3 hierarchy selects the *filtered* gene models only - JGI ``jat_label``
``genes_filtered`` (i.e. the GeneCatalog / FilteredModels ``.gff``) - never the
unfiltered ``genes_all`` models, exactly as the XML parser did.

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
import subprocess
import requests
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from urllib.parse import unquote
from mycotools.lib.kontools import format_path, outro, intro, setup_logging
from mycotools.lib.dbtools import loginCheck
from pathlib import Path

logger = logging.getLogger(__name__)

# JGI Data Portal API endpoints (non-Globus)
SIGNON_URL = "https://signon.jgi.doe.gov/signon/create"
SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
RESTORE_URL = "https://files.jgi.doe.gov/request_archived_files/"
DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/"


def jgi_login(user, pwd):
    """Login via JGI's prescribed method by creating a cookie cache and
    downloading JGI's sign-in file."""

    null = str(Path("~/.nulljgi_dwnld").expanduser())

    login_cmd = subprocess.call(
        [
            "curl",
            "https://signon.jgi.doe.gov/signon/create",
            "--data-urlencode",
            "login=" + str(user),
            "--data-urlencode",
            "password=" + str(pwd),
            "-c",
            "cookies",
            "-o",
            null,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    return login_cmd


def dwnld_xml(output, ome, max_tempts=2):
    attempts = 0
    while not Path(f"{output}/{ome}.xml").is_file() and attempts < max_tempts:
        attempts += 1
        xml_cmd = subprocess.call(
            [
                "curl",
                "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism="
                + str(ome),
                "-b",
                "cookies",
                "-o",
                f"{output}/{ome}.xml",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if xml_cmd != 0:
            logger.error(f"\t{ome} xml curl error: {xml_cmd}")
    if not Path(f"{output}/{ome}.xml").is_file():
        return -1
    else:
        return xml_cmd


def is_directory_xml(xml_data):
    """Return whether `xml_data` is a well-formed JGI organism-directory XML
    document rather than an HTML error/redirect page.

    JGI has, at times, served 302 redirect pages or error HTML (e.g. when an
    endpoint is deprecated) in place of the directory XML. Those must be caught
    here so they never reach the XML parser, which would otherwise raise an
    unhandled ParseError and abort the entire download run."""
    if not xml_data or not xml_data.strip():
        return False
    head = xml_data.lstrip().lower()
    if head.startswith("<!doctype html") or head.startswith("<html"):
        return False
    try:
        ET.fromstring(xml_data)
    except ET.ParseError:
        return False
    return True


def retrieve_xml(ome, output):
    """Retrieve JGI xml file tree. First check if it already exists, if not then
    download it using JGI's prescribed method. Then open the xml and check for
    the common 'Portal does not exist' error, and validate that the response is
    actually XML (not an HTML error/redirect page). Report failures so the ome
    is skipped rather than crashing the parser."""

    xml_path = f"{output}/{ome}.xml"

    if not Path(xml_path).exists():
        dwnld_xml(output, ome)
    if not Path(xml_path).exists():
        return None  # curl produced no file; caller will retry

    with open(xml_path, "r") as xml_raw:
        xml_data = xml_raw.read()

    if xml_data == "Portal does not exist":
        logger.error("\t`" + ome + " not in JGIs `organism` database")
        Path(xml_path).unlink()
        return 1
    if not xml_data:
        Path(xml_path).unlink()
        return None
    if not is_directory_xml(xml_data):
        # non-XML response (e.g. an HTML error/redirect page from a deprecated
        # JGI endpoint); discard so it never reaches parse_xml
        logger.error(
            f"\t`{ome}` returned a non-XML directory response "
            "(the JGI download API may have changed); skipping"
        )
        Path(xml_path).unlink()
        return 1

    return -1


def parse_xml(ft, xml_file, masked=False, forbidden={}, filtered=True):
    """Parse the XML data to obtain the file types of interest based on
    predefined hashes that contain the known subdirectories associated with JGI
    organization"""

    if ft == "fna":
        if masked:
            ft += "$masked"
        else:
            ft += "$unmasked"

    # set the initial hashes for the XML hierarchy - relate file types to their
    # hierarchy structure
    ft2xt = {
        "fna$masked": {"assembly"},
        "fna$unmasked": {"assembly"},
        "gff": {"annotation"},
        "gff3": {"annotation"},
        "transcripts": {"annotation"},
        "est": {"ests and est clusters", "transcriptome"},
    }
    ft2fh = {
        "fna$masked": ["genome assembly (masked)", "assembled scaffolds (masked)"],
        "fna$unmasked": [
            "assembled scaffolds (unmasked)",
            "genome assembly (unmasked)",
        ],
        "gff": ["genes"],
        "gff3": ["genes"],
        "transcripts": ["transcripts"],
        "est": ["ests", "transcriptome assembly"],
    }

    ft2fn = {
        "fna$masked": ["masked", "Genome Assembly (masked)"],
        "fna$unmasked": [
            "AssembledScaffolds",
            "scaffolds",
            "Genome Assembly (unmasked)",
            "AssemblyScaffolds",
        ],
        "gff3": ["GeneCatalog", "FilteredModels"],
        "transcripts": ["transcripts"],
        "est": ["EST"],
    }
    ft2fe = {
        "fna$masked": {"fasta", "fa", "fna", "fsa"},
        "fna$unmasked": {"fasta", "fa", "fna", "fsa"},
        "gff": {"gff", "gff3"},
        "gff3": {"gff", "gff3"},
        "transcripts": {"fa", "fasta", "fna", "fsa"},
        "est": {"fa", "fasta", "fna", "fsa"},
    }

    url, md5, filename = None, False, None

    # parse the XML file; a malformed/non-XML file (e.g. an HTML error page that
    # slipped through) must not abort the whole run
    try:
        tree = ET.parse(xml_file)
    except ET.ParseError as parse_error:
        logger.error(f"\tmalformed JGI directory XML {xml_file}: {parse_error}")
        return None, None, False, None
    root = tree.getroot()
    flip = True
    org_name = None
    has_flipped = False
    attempt = 0

    # flip is a way to rerun the loop if the file type changes (e.g. from
    # masked to unmasked); parse through the XML hiearchy in accord with the
    # hashes established above
    while flip and attempt < 10:
        attempt += 1
        for child in root:
            # conserved subdirectory we need
            if "Files" == child.attrib["name"]:
                for chil1 in child:
                    # does this subdirectory match what we need for our
                    # filetype?
                    if chil1.attrib["name"].lower() in ft2xt[ft]:
                        # we only want the filtered models for annotations/RNA
                        if ft in {"gff3", "gff", "transcripts", "est"}:
                            for chil2 in chil1:
                                if (
                                    chil2.attrib["name"]
                                    .lower()
                                    .startswith("filtered models (")
                                ):
                                    chil1 = chil2
                                    break
                        # continue parsing toward the files of interest
                        for chil2 in chil1:
                            if any(
                                x == chil2.attrib["name"].lower() for x in ft2fh[ft]
                            ):
                                for chil3 in chil2:
                                    t_url = chil3.attrib["url"]
                                    try:
                                        org_name = chil3.attrib["label"]
                                    except KeyError:
                                        pass
                                    # we want to avoid tape files as we cannot
                                    # download them readily
                                    if (
                                        "get_tape_file" not in t_url
                                        and t_url not in forbidden
                                    ):
                                        if all(
                                            x not in chil3.attrib["filename"]
                                            for x in ft2fn[ft]
                                        ):
                                            continue
                                        file_ext_srch = re.search(
                                            r"\.([^\.]+)$", chil3.attrib["filename"]
                                        )
                                        if file_ext_srch is not None:
                                            file_ext = file_ext_srch[1]
                                            if file_ext == "gz":
                                                file_ext_srch = re.search(
                                                    r"\.([^\.]+)\.gz$",
                                                    chil3.attrib["filename"],
                                                )
                                                if file_ext_srch is not None:
                                                    file_ext = file_ext_srch[1]
                                            if file_ext not in ft2fe[ft]:
                                                continue
                                        else:
                                            continue

                                        url = chil3.attrib["url"]
                                        filename = chil3.attrib["filename"]
                                        # all requirements satisfied
                                        try:
                                            md5 = chil3.attrib["md5"]
                                            break
                                        # continue on to find an md5, or omit
                                        # if exhaustively searched
                                        except KeyError:
                                            pass

        # if the unmasked genome is not present, then query for the masked and
        # vice versa
        if not url and ft == "fna$masked":
            ft = "fna$unmasked"
            if not has_flipped:
                flip = True
            else:
                flip = False
        elif not url and ft == "fna$unmasked":
            ft = "fna$masked"
            if not has_flipped:
                flip = True
            else:
                flip = False
        else:
            break

    return filename, url, md5, org_name


def handle_redirect_307(
    dwnld_data, dwnld, dwnld_url, file_type, xml_file, masked, url, urls, spacer
):
    """Handle a redirection error by identifying a new file URL to download
    from, or return the original if none exist"""
    logger.info(
        spacer + "\t" + dwnld + " link has moved. " + "Trying a different link."
    )
    filename, n_url, dwnld_md5, t_org_name = parse_xml(
        file_type, xml_file, masked=masked, forbidden={url}.union(urls)
    )

    if n_url:
        url = n_url
        dwnld_url = prefix + url.replace("&amp;", "&")
        dwnld = f"{output}{file_type}/{Path(dwnld_url).name}"
    return url, dwnld_url, dwnld, {url}.union(urls), t_org_name


def no_md5_checks(dwnld, md5, spacer):
    """If there is no MD5, simply check the file has content in it"""
    check_size = subprocess.run(["wc", "-l", dwnld], stdout=subprocess.PIPE)
    check_size_res = check_size.stdout.decode("utf-8")
    logger.info(spacer + "\t\tFile exists - no md5 to check.")
    check_size_find = re.search(r"\d+", check_size_res)
    size = check_size_find[0]
    if int(size) < 10:
        logger.warning(spacer + "\tInvalid file size.")
    else:
        md5 = None
    return md5


def jgi_dwnld(ome, file_type, output, masked=True, spacer="\t"):
    """Download JGI files. For each type of file, use regular expressions to
    gather the URL from the file. Grab the md5checksum if possible from the
    xml as well. Create arbitrary values for md5 and curl_cmd. If the download file
    already exists, then run an md5 checksum if an md5 value exists in the xml.

    If the md5 does not match then open and check for typical errors. If those
    errors exist, obtain the alternative download URL from the xml. If the md5
    matches, pass through the rest of the function.

    If there is no download md5, then check to see if the file is greater than 10
    lines as a proxy to make sure that the file isn't empty/blatantly wrong. If it
    passes this test, change the values to not enter the while loop at the end of
    the function.

    The while loop following the file exists error allows for 3 attempts. If a file
    is downloaded it will check its md5 using the xml reference. If it fails, it will
    proceed via the error checking above, wait a minute, and reattempt. If there is no
    download md5 from the xml it will proceed via the error checking above as well. If it
    passes the line count check, it will exit the loop - otherwise, it will attempt to
    gather a new URL and restart after another minute wait. This is an unfortunate
    circumnavigation of JGI's cryptic maximum ping / time before booting."""

    # prepare data structures
    prefix = "https://genome.jgi.doe.gov"
    xml_file = f"{output}xml/{ome}.xml"
    preexisting, check = False, 1

    # stop gap for legacy input
    if file_type == "gff":
        file_type += "3"

    ran_dwnld = False
    org_name = None

    # acquire the filename, URL, and MD5 from the xml for the file type of
    # interest
    filename, url, dwnld_md5, t_org_name = parse_xml(file_type, xml_file, masked=masked)
    if t_org_name:
        org_name = t_org_name
    if not dwnld_md5:
        dwnld_md5 = None

    # if there is a URL present, begin the downloading process
    if url:
        f_urls = {url}
        md5 = False
        attempt = 0
        curl_cmd = 420

        dwnld_url = prefix + url.replace("&amp;", "&")

        dwnld = f"{output}{file_type}/{Path(dwnld_url).name}"
        unzip_dwnld = re.sub(r"\.gz$", "", dwnld)
        # assume unzipped downloads have passed the checks
        if Path(unzip_dwnld).is_file():
            md5 = dwnld_md5
            curl_cmd = 0
            check = unzip_dwnld
            preexisting = True

        # if the file currently exists, then check its MD5
        elif Path(dwnld).exists():
            if dwnld_md5:
                md5_cmd = subprocess.run(
                    ["md5sum", dwnld], stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                md5_res = md5_cmd.stdout.decode("utf-8")
                md5_find = re.search(r"\w+", md5_res)
                md5 = md5_find[0]

            if md5 == dwnld_md5:
                curl_cmd = 0
                check = dwnld
                preexisting = True
            # if the MD5 does not equal the download MD5 then check the file
            else:
                while True:
                    try:
                        with open(dwnld, "r") as dwnld_data_raw:
                            dwnld_data = dwnld_data_raw.read()
                        if re.search("307 Temporary Redirect", dwnld_data):
                            url, dwnld_url, dwnld, f_urls, t_org_name = (
                                handle_redirect_307(
                                    dwnld_data,
                                    dwnld,
                                    dwnld_url,
                                    file_type,
                                    xml_file,
                                    masked,
                                    url,
                                    f_urls,
                                    spacer,
                                )
                            )
                            if t_org_name:
                                org_name = t_org_name
                        break
                    except FileNotFoundError:
                        md5 = False
                        break
                    except UnicodeDecodeError:
                        if not dwnld_md5:
                            md5 = no_md5_checks(dwnld, md5, spacer)
                            if not md5:
                                preexisting = True
                                check = dwnld
                        else:
                            logger.warning(spacer + "\tmd5 does not match.")
                        break

        # while the MD5 doesn't match, or there is a curl error, try up to 3
        # times to download the file
        while md5 != dwnld_md5 and curl_cmd != 0 and attempt < 3:
            attempt += 1
            curl_cmd = subprocess.call(
                ["curl", dwnld_url, "-b", "cookies", "-o", dwnld],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            ran_dwnld = True

            if curl_cmd == 0:
                check = dwnld

                # acquire the MD5
                if dwnld_md5:
                    md5_cmd = subprocess.run(
                        ["md5sum", dwnld],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    )
                    md5_res = md5_cmd.stdout.decode("utf-8")
                    md5_find = re.search(r"\w+", md5_res)
                    try:
                        md5 = md5_find[0]
                    except TypeError:
                        md5 = False
                        if not Path(dwnld).is_file():
                            attempt += 1
                            continue
                # if there is no MD5 attempt the crude file check
                if not dwnld_md5:
                    try:
                        with open(dwnld, "r") as dwnld_data_raw:
                            dwnld_data = dwnld_data_raw.read()
                        if re.search("307 Temporary Redirect", dwnld_data):
                            t_url, dwnld_url, dwnld, f_urls, t_org_name = (
                                handle_redirect_307(
                                    dwnld_data,
                                    dwnld,
                                    dwnld_url,
                                    file_type,
                                    xml_file,
                                    masked,
                                    url,
                                    f_urls,
                                    spacer,
                                )
                            )
                            if t_org_name:
                                org_name = t_org_name

                            if t_url == url:
                                logger.warning(spacer + "\t\tNo valid alternative")
                                attempt = 4
                                break
                            else:
                                url = t_url
                        else:
                            md5 = no_md5_checks(dwnld, md5, spacer)
                            if not md5:
                                break
                    except FileNotFoundError:
                        pass
                    except UnicodeDecodeError:
                        pass
                    # this is slow, and was arbitrarily set to not overping JGI
                    time.sleep(60)

                # the download may have failed, so prepare to retry
                elif md5 != dwnld_md5 and attempt == 1:
                    logger.error(
                        f"{spacer}\tmd5 does not match JGI. " + f"Attempt {attempt}"
                    )
                    curl_cmd = -1
                    check = 2
                    while True:
                        try:
                            with open(dwnld, "r") as dwnld_data_raw:
                                dwnld_data = dwnld_data_raw.read()
                            if re.search("307 Temporary Redirect", dwnld_data):
                                t_url, dwnld_url, dwnld, f_urls, t_org_name = (
                                    handle_redirect_307(
                                        dwnld_data,
                                        dwnld,
                                        dwnld_url,
                                        file_type,
                                        xml_file,
                                        masked,
                                        url,
                                        f_urls,
                                        spacer,
                                    )
                                )
                                if t_org_name:
                                    org_name = t_org_name
                                if t_url == url:
                                    logger.warning(spacer + "\t\tNo valid alternative")
                                    attempt = 4
                                    break
                            break
                        except FileNotFoundError:
                            pass
                            break
                        except UnicodeDecodeError:
                            pass
                            break
                        time.sleep(60)
                # if there are two fails, attempt a new URL
                elif md5 != dwnld_md5 and attempt == 2:
                    logger.error(f"{spacer}\tmd5 does not match JGI. Attempt {attempt}")
                    curl_cmd = -1
                    filename, n_url, dwnld_md5, t_org_name = parse_xml(
                        file_type, xml_file, masked=masked, forbidden=f_urls
                    )
                    if t_org_name:
                        org_name = t_org_name

                    if n_url:
                        url = n_url
                        dwnld_url = prefix + url.replace("&amp;", "&")
                        f_ulrs = {url}.union(f_urls)
                        dwnld = f"{output}{file_type}/{Path(dwnld_url).name}"
                    time.sleep(60)
                    check = 2
                elif md5 != dwnld_md5:
                    logger.error(f"{spacer}\tmd5 does not match JGI. Attempt {attempt}")
                    check = 2
            else:
                logger.error(
                    f"{spacer}\tFailed to retrieve {file_type}. `curl` error: "
                    + f"{curl_cmd}\n{spacer}\tAttempt {attempt}"
                )
                check = 2

        # three strikes and the file is out
        if attempt == 3:
            if md5 != dwnld_md5:
                logger.warning(spacer + "\tExcluding from database - potential failure")
                curl_cmd = 0
            if curl_cmd != 0:
                logger.error(spacer + "\tFile failed to download")

    return check, preexisting, file_type, ran_dwnld, org_name


# ===========================================================================
# JGI Data Portal API - non-Globus download implementation (see module docstring)
# ===========================================================================

# Ordered, most-preferred-first jat_labels and acceptable (un-gzipped) file
# formats per download type. This encodes the same hierarchy as parse_xml's
# ft2xt/ft2fh/ft2fn/ft2fe hashes above - most importantly, "gff3" resolves only
# to the filtered gene models (genes_filtered), never genes_all.
_TYPE_LABELS = {
    "gff3": (["genes_filtered"], {"gff", "gff3"}),
    "faa": (["proteins_filtered"], {"fasta", "fa", "aa"}),
    "transcript": (["transcripts_filtered"], {"fasta", "fa", "fna", "fsa", "nt"}),
    "est": (["ests", "est_clusters"], {"fasta", "fa", "fna", "fsa"}),
}


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


def select_file(files, ftype, masked=True):
    """Choose the single best file record for `ftype` from an organism's file
    list, mirroring parse_xml's selection hierarchy.

    Preference, in order:
      1. immediately-available (RESTORED) files over archived (PURGED) ones -
         matching the legacy parser's avoidance of on-tape ``get_tape_file``
         URLs (and its masked->unmasked flip when the preferred assembly was on
         tape);
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
    (genus, species, strain), dropping a trailing version token. Mirrors the
    legacy label parse (strain is the remaining words concatenated)."""
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
    restore, and download endpoints. Exits (100) after repeated failures, as the
    legacy login did."""
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
                logger.warning(
                    f"{spacer}\t{portal_id} search error (attempt {attempt}): {error}"
                )
                time.sleep(2)
                continue
            if resp.status_code != 200:
                logger.warning(f"{spacer}\t{portal_id} search HTTP {resp.status_code}")
                time.sleep(2)
                continue
            try:
                data = resp.json()
            except ValueError:
                logger.warning(f"{spacer}\t{portal_id} search returned non-JSON")
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
    the timeout / expiry is reached (returns False).

    Progress is logged so the user can see the tape->disk transfer advance:
    every status transition (queued -> retrieving -> staging -> ready) is
    reported, plus a heartbeat every `heartbeat` seconds while a stage lingers."""
    if not status_url:
        return False
    tag = f"{label}: " if label else ""
    waited = 0
    last_status = None
    last_heartbeat = 0
    while waited <= timeout:
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
    restore_timeout=60,
    poll_interval=30,
    request_delay=3,
):
    """Download MycoCosm data for the JGI portal ids in `df` via the JGI Data
    Portal API (non-Globus), preserving the legacy contract: `df` gains
    ``<type>_path`` columns (e.g. fna_path, gff3_path) plus genus/species/strain,
    and the function returns (df, failed_portal_ids).

    Archived (PURGED) files are restored from tape before download; if a restore
    does not finish within `restore_timeout` the portal id is deferred (added to
    the returned failure set) so a later rerun can pick it up once ready."""
    if "assembly_acc" in df.columns:
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
    for i, row in tqdm(df.iterrows(), total=len(df)):
        portal_id = row[ome_col]

        org, files = search_organism(session, portal_id, spacer=spacer)
        if not org:
            logger.warning(f"{spacer}\t{portal_id} not found in JGI MycoCosm")
            ome_set.add(portal_id)
            continue

        org_id = org.get("id")
        top_hit = (org.get("top_hit") or {}).get("_id")
        portal = org.get("mycocosm_portal_id") or portal_id

        selected = {}
        for typ in dwnlds:
            chosen = select_file(files, typ, masked=masked)
            if chosen is not None:
                selected[typ] = chosen
        if not selected:
            logger.warning(f"{spacer}\t{portal_id}: no target files available")
            ome_set.add(portal_id)
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
            logger.warning(
                f"{spacer}\t{portal_id}: {len(on_tape)} file(s) are archived on TAPE and "
                "must be restored to disk before download - "
                + ", ".join(f.get("file_name", f["_id"]) for f in on_tape)
            )
            status_url = request_restore(
                session,
                token,
                _mycocosm_ids(org_id, top_hit, portal, [f["_id"] for f in on_tape]),
                spacer=spacer,
            )
            if not poll_restore(
                session,
                status_url,
                timeout=restore_timeout,
                interval=poll_interval,
                spacer=spacer,
                label=portal_id,
            ):
                logger.warning(
                    f"{spacer}\t{portal_id}: tape restore still pending; deferring this "
                    "genome (rerun later to resume once JGI has staged it to disk)"
                )
                ome_set.add(portal_id)
                continue

        file_ids = [f["_id"] for f in selected.values()]
        dest_zip = os.path.join(tmp_dir, f"{portal_id}.zip")
        if not download_zip(
            session,
            token,
            _mycocosm_ids(org_id, top_hit, portal, file_ids),
            dest_zip,
            spacer=spacer,
        ):
            logger.warning(f"{spacer}\t{portal_id}: download failed")
            ome_set.add(portal_id)
            continue

        wanted, type_dest = {}, {}
        for typ, f in selected.items():
            name = f["file_name"]
            dest = os.path.join(output, typ, name)
            wanted[name] = dest
            type_dest[typ] = (name, dest)
        extracted = extract_zip(dest_zip, wanted, spacer=spacer)
        if Path(dest_zip).is_file():
            Path(dest_zip).unlink()

        essential_failed = False
        for typ, (name, dest) in type_dest.items():
            if name in extracted and Path(dest).is_file():
                df.at[i, typ + "_path"] = dest
                logger.info(f"{spacer}\t{portal_id} {typ}: {name}")
            else:
                logger.warning(f"{spacer}\t{portal_id}: {typ} missing from archive")
                if typ in ("fna", "gff3"):
                    essential_failed = True

        genus, species, strain = parse_org_name(org.get("name") or "")
        _fill_if_empty(df, i, "genus", genus)
        _fill_if_empty(df, i, "species", species)
        _fill_if_empty(df, i, "strain", strain)

        if essential_failed:
            ome_set.add(portal_id)
        if request_delay:
            time.sleep(request_delay)

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
    parser.add_argument("-o", "--output", default=str(Path.cwd()), help="Output dir")
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

    ncbi_email, ncbi_api, user, pwd = loginCheck(ncbi=False)

    args_dict = {
        "JGI Table": args.input,
        "Assemblies": args.assembly,
        "RepeatMasked": not args.nonmasked,
        "Proteomes": args.proteome,
        ".gff3's": args.gff,
        "Transcripts": args.transcript,
        "EST": args.est,
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
    )
    jgi_df = jgi_df.rename(columns={"assembly_acc": "#assembly_acc"})
    jgi_df["source"] = "jgi"
    jgi_df["restriction"] = "no"
    jgi_df.to_csv(str(Path(args.input)) + ".predb.tsv", sep="\t", index=False)

    outro(start_time)


if __name__ == "__main__":
    cli()

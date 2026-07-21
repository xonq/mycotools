#! /usr/bin/env python3

# NEED to make login check more intuitive and easier for NCBI only

from __future__ import annotations

import os
import re
import sys
import copy
import json
import logging
import time
import base64
import urllib
import getpass
import zipfile
import datetime
import subprocess
import random
from tqdm import tqdm
from Bio import Entrez
from io import StringIO
from typing import Any, Dict, Iterable, Mapping, Optional, Union
from collections import defaultdict
from mycotools.lib.kontools import (
    collect_files,
    format_path,
    read_json,
    write_json,
)
from pathlib import Path

logger = logging.getLogger(__name__)


class mtdb(dict):
    """
    MycotoolsDB (mtdb) class. Designed to support high throughput "pandas-like"
    interface without the overhead of pandas.
    """

    columns = [
        "ome",
        "genus",
        "species",
        "strain",
        "taxonomy",
        "version",
        "source",
        "biosample",
        "assembly_acc",
        "acquisition_date",
        "published",
        "fna",
        "faa",
        "gff3",
    ]

    #: columns that may serve as the MTDB index (unique, hashable keys)
    _index_columns = frozenset({"assembly_acc", "ome"})
    #: taxonomic ranks excluded when assimilating taxonomy dictionaries
    _forbidden_tax_ranks = frozenset(
        {
            "no rank",
            "subkingdom",
            "genus",
            "species",
            "species group",
            "varietas",
            "forma",
        }
    )

    # NEED a detect index feature for adding dicts in with alternative indices
    def __init__(
        self,
        db: Union[str, Mapping[str, Any], "mtdb", None] = None,
        index: Optional[str] = None,
        add_paths: bool = True,
    ):
        if not db:
            if not index:
                super().__init__({x: [] for x in self.columns})
            else:
                super().__init__()
        elif isinstance(db, dict):
            if not index:
                if all(isinstance(db[x], list) for x in self.columns):
                    super().__init__(db)
                    self.index = None
            else:
                try:
                    if set(db[list(db.keys())[0]].keys()).union({index}) == set(
                        self.columns
                    ):
                        super().__init__(db)
                except AttributeError:
                    if set(db[list(db.keys())[0]][0]).union({index}) == set(
                        self.columns
                    ):
                        super().__init__(db)
        else:
            super().__init__(self.db2df(db, add_paths=add_paths))
        self.index = index

    def mtdb2pd(self) -> "pd.DataFrame":
        import pandas as pd

        copy_mtdb = copy.deepcopy(self)
        copy_mtdb = copy_mtdb.reset_index()
        return pd.DataFrame(copy_mtdb)  # assume pd is imported

    @staticmethod
    def pd2mtdb(df: "pd.DataFrame") -> "mtdb":  # legacy integration
        df = df.fillna("")
        for i, row in df.iterrows():
            if not row["gff3"]:
                row["gff3"] = os.environ["MYCOGFF3"] + row["ome"] + ".gff3"
            if not row["faa"]:
                row["faa"] = os.environ["MYCOFAA"] + row["ome"] + ".faa"
            if not row["fna"]:
                row["fna"] = os.environ["MYCOFNA"] + row["ome"] + ".fna"
        db = mtdb({x: list(df[x]) for x in mtdb.columns})
        return db

    def db2df(self, db_path: str, add_paths: bool = True) -> Dict[str, list]:
        df = defaultdict(list)
        if Path(db_path).stat().st_size == 0:
            return {x: [] for x in mtdb.columns}
        with open(format_path(db_path), "r") as raw:
            data = [
                x.rstrip().split("\t")
                for x in raw
                if not x.startswith("#") and x.rstrip()
            ]
        columns = self.columns
        for entry in data:
            [df[c].append("") for c in columns]  # add a blank entry to each
            # column
            for i, d in enumerate(entry):
                df[columns[i]][-1] = d
            try:
                df["taxonomy"][-1] = self.read_tax(df["taxonomy"][-1])
            except json.decoder.JSONDecodeError:
                logger.error("malformed taxonomy: %s", df["taxonomy"][-1])
                sys.exit()
            df["taxonomy"][-1]["genus"] = df["genus"][-1]
            df["taxonomy"][-1]["species"] = df["genus"][-1] + " " + df["species"][-1]
            df["taxonomy"][-1]["strain"] = df["strain"][-1]

        if not add_paths:
            return df
        try:
            for i, ome in enumerate(df["ome"]):
                if not df["fna"][i]:
                    if {"MYCOFNA", "MYCOFAA", "MYCOGFF3"}.difference(
                        set(os.environ.keys())
                    ):
                        raise FileNotFoundError(
                            "You are not connected to a primary MTDB. "
                            + "Standalone databases need absolute paths"
                        )
                    df["fna"][i] = os.environ["MYCOFNA"] + ome + ".fna"
                    df["faa"][i] = os.environ["MYCOFAA"] + ome + ".faa"
                    df["gff3"][i] = os.environ["MYCOGFF3"] + ome + ".gff3"
                elif df["fna"][i] == ome + ".fna":
                    if {"MYCOFNA", "MYCOFAA", "MYCOGFF3"}.difference(
                        set(os.environ.keys())
                    ):
                        raise FileNotFoundError(
                            "You are not connected to a primary MTDB. "
                            + "Standalone databases need absolute paths"
                        )
                    df["fna"][i] = os.environ["MYCOFNA"] + ome + ".fna"
                    df["faa"][i] = os.environ["MYCOFAA"] + ome + ".faa"
                    df["gff3"][i] = os.environ["MYCOGFF3"] + ome + ".gff3"
        except KeyError:
            logger.error("MycotoolsDB not in path, cannot delineate biofile paths")

        return df

    def df2db(
        self,
        db_path: Optional[str] = None,
        headers: bool = False,
        paths: bool = False,
    ) -> None:
        df = copy.copy(self)
        df = df.reset_index()
        output = mtdb(
            {
                k: v
                for k, v in sorted(self.set_index("ome").items(), key=lambda x: x[0])
            },
            index="ome",
        )
        # does this work if its not an inplace change
        abb_paths = {
            "faa": [os.environ["MYCOFAA"], ".faa"],
            "fna": [os.environ["MYCOFNA"], ".fna"],
            "gff3": [os.environ["MYCOGFF3"], ".gff3"],
        }
        if db_path:
            with open(db_path, "w") as out:
                if headers:
                    out.write("#" + "\t".join(self.columns) + "\n")
                for ome in output:
                    if not paths:
                        for file_type in ["fna", "faa", "gff3"]:
                            output[ome][file_type] = output[ome][file_type].replace(
                                abb_paths[file_type][0] + ome + abb_paths[file_type][1],
                                "",
                            )  # abbreviate when possible
                    for rank in ["species", "genus", "strain"]:
                        try:
                            del output[ome]["taxonomy"][rank]
                        except (KeyError, TypeError) as e:
                            pass

                    if output[ome]["taxonomy"]:
                        output[ome]["taxonomy"] = json.dumps(output[ome]["taxonomy"])
                    else:
                        output[ome]["taxonomy"] = "{}"
                    if not output[ome]["published"]:
                        output[ome]["published"] = ""
                    out.write(
                        ome
                        + "\t"
                        + "\t".join([str(output[ome][x]) for x in output[ome]])
                        + "\n"
                    )
        else:
            if headers:
                print("#" + "\t".join(self.columns), flush=True)

            for ome in output:
                if not paths:
                    for file_type in ["fna", "faa", "gff3"]:
                        output[ome][file_type] = output[ome][file_type].replace(
                            abb_paths[file_type][0] + ome + abb_paths[file_type][1], ""
                        )  # abbreviate when possible
                for rank in ["species", "genus", "strain"]:
                    try:
                        del output[ome]["taxonomy"][rank]
                    except (KeyError, TypeError) as e:
                        pass
                output[ome]["taxonomy"] = json.dumps(output[ome]["taxonomy"])
                print(
                    ome + "\t" + "\t".join([str(output[ome][x]) for x in output[ome]]),
                    flush=True,
                )

    def set_index(self, column: Optional[str] = "ome", inplace: bool = False) -> "mtdb":
        data, retry, error, df, columns = (
            {},
            bool(column),
            False,
            copy.copy(self),
            copy.copy(self.columns),
        )
        if not column:
            return df.reset_index()
        elif column not in self._index_columns:
            raise KeyError(f'MTDB index must be "assembly_acc"/"ome"')
        elif df.index and not df.keys():  # empty df
            return mtdb({}, index=column)
        elif not df.index:
            if not df["ome"]:
                return mtdb({}, index=column)
        while retry:
            oldCol = set()
            try:
                columns.pop(columns.index(column))
            except (ValueError, IndexError) as e:
                df = df.reset_index()  # will this actually reset the index
                columns.pop(columns.index(column))
            try:
                df[column]
            except KeyError:
                df = df.reset_index()
            if len(df[column]) == 0:
                return mtdb({}, index=column)
            if column in {
                "assembly_acc",
                "ome",
                "fna",
                "gff3",
                "faa",
            }:  # if it's unique
                for i, v in enumerate(df[column]):
                    data[v] = {}
                    for head in columns:
                        data[v][head] = df[head][i]
            else:
                for i, v in enumerate(df[column]):
                    if v not in data:
                        data[v] = []
                    data[v].append({})
                    for head in columns:
                        try:
                            data[v][-1][head] = df[head][i]
                        except KeyError:  # if an index exists
                            if error:
                                logger.error("invalid column")
                                return self
                            df = df.reset_index()
                            error = True
                            break

            retry = False
        return mtdb(data, column)

    def reset_index(self) -> "mtdb":
        df = copy.copy(self)
        if df.index:
            data = {x: [] for x in mtdb().columns}
            data[df.index] = list(df.keys())
            if df.index in {"assembly_acc", "ome", "biosample", "fna", "gff3", "faa"}:
                for key in df:
                    for otherCol in df[key]:
                        data[otherCol].append(df[key][otherCol])
            else:
                for key in df:
                    for v in df[key]:
                        for otherCol in v:
                            data[otherCol].append(v[otherCol])
            return mtdb(data, index=None)
        else:
            return df

    def append(self, info: Optional[Mapping[str, Any]] = None) -> "mtdb":
        #        if any(x not in set(self.columns) for x in info):
        #           raise KeyError('Invalid keys: ' + str(set(info.keys()).difference(set(self.columns))))
        if info is None:
            info = {}
        index = self.index
        df = copy.copy(self)
        df = df.reset_index()
        info = {
            **info,
            **{k: None for k in set(self.columns).difference(set(info.keys()))},
        }
        for key in self.columns:
            df[key].append(info[key])
        return df.set_index(index)

    @staticmethod
    def read_tax(
        taxonomy_string: Union[str, Mapping[str, Any], None]
    ) -> Dict[str, Any]:
        """Read taxonomy from an MTDB by converting the string into a dictionary"""
        tax_strs = [
            "superkingdom",
            "kingdom",
            "phylum",
            "subphylum",
            "class",
            "order",
            "family",
            "subfamily",
        ]
        if taxonomy_string:
            if isinstance(taxonomy_string, str):
                dict_string = taxonomy_string.replace("'", '"')
                try:
                    tax_dict = json.loads(dict_string)
                except TypeError:
                    tax_dict = {}
            else:
                tax_dict = taxonomy_string
            try:
                tax_dict = {
                    **tax_dict,
                    **{x: "" for x in tax_strs if x not in tax_dict},
                }
            except TypeError:  # inappropriate tax_dict in the column
                tax_dict = {x: "" for x in tax_strs}
            return tax_dict
        else:
            return {}

    @staticmethod
    def _reconcile_tax_dicts(
        genera: Iterable[str],
        tax_dicts: Mapping[str, Mapping[str, Any]],
        forbid: Iterable[str],
    ) -> Dict[str, Dict[str, Any]]:
        """Drop empty/forbidden taxonomy entries and backfill missing genera"""
        forbid = set(forbid)
        tax_dicts = {x: tax_dicts[x] for x in tax_dicts if tax_dicts[x]}
        for genus in tax_dicts:
            tax_dicts[genus] = {
                rank: name
                for rank, name in tax_dicts[genus].items()
                if rank not in forbid
            }
        for miss in set(genera).difference(tax_dicts.keys()):
            tax_dicts[miss] = {}
        return tax_dicts

    def prepare_tax_dicts(
        self, tax_dicts: Optional[Dict[str, Any]] = None
    ) -> "tuple[set, Dict[str, Any]]":
        """Identify the genera that do not have higher taxonomy ascribed to them"""
        if tax_dicts is None:
            tax_dicts = {}
        need_tax = set()
        for ome, entry in self.set_index("ome").items():
            if entry["genus"] in tax_dicts:
                continue
            tax_json = self.read_tax(entry["taxonomy"])
            if any(
                name
                for rank, name in tax_json.items()
                if rank not in {"genus", "species", "strain"}
            ):
                tax_dicts[entry["genus"]] = tax_json
            else:
                need_tax.add(entry["genus"])
        need_tax = need_tax.difference(tax_dicts.keys())
        return need_tax, tax_dicts

    def assimilate_tax(
        self,
        tax_dicts: Mapping[str, Mapping[str, Any]],
        forbid: Optional[Iterable[str]] = None,
    ) -> "tuple[mtdb, Dict[str, Any]]":
        """Assign the resolved taxonomy dictionaries to this MTDB's taxonomy column"""
        if forbid is None:
            forbid = self._forbidden_tax_ranks
        tax_dicts = self._reconcile_tax_dicts(set(self["genus"]), tax_dicts, forbid)
        for i, genus in enumerate(self["genus"]):
            self["taxonomy"][i] = tax_dicts[genus]
        return mtdb(self), tax_dicts

    def infer_rank(self, lineage: str) -> str:
        """Identify the taxonomic rank associated with an inputted lineage of
        interest"""
        linlow, rank = lineage.lower(), None
        for ome, row in self.items():
            if linlow in set([x.lower() for x in row["taxonomy"].values()]):
                rev_dict = {
                    k.lower(): v
                    for k, v in zip(row["taxonomy"].values(), row["taxonomy"].keys())
                }
                rank = rev_dict[lineage]

        if not rank:
            raise KeyError(f"no entry for {lineage}")

        return rank

    def extract_unique(self, allowed: int = 1, rank: str = "species") -> "mtdb":
        """Extract unique rank from an MTDB"""
        keys = copy.deepcopy(list(self.keys()))
        random.shuffle(keys)
        prep_db0 = {x: self[x] for x in keys}
        prep_db1 = mtdb().set_index("ome")
        if rank == "strain":
            found = set()
            for ome, row in prep_db0.items():
                name = row["taxonomy"]["species"] + " " + row["strain"]
                if name not in found:
                    prep_db1[ome] = row
                found_prep = list(found)
                found_prep.append(name)
                found = set(found_prep)
        else:
            found = defaultdict(int)
            for ome, row in prep_db0.items():
                name = row["taxonomy"][rank]
                found[name] += 1
                if found[name] <= allowed:
                    prep_db1[ome] = row

        return prep_db1

    def extract_tax(self, lineages: Union[str, Iterable[str]]) -> "mtdb":
        """Extract specific taxonomic lineages of interest based on their rank"""
        if isinstance(lineages, str):
            lineages = [lineages]
        lineages = set(x.lower() for x in lineages)
        rank_dict = {k: self.infer_rank(k) for k in list(lineages)}
        ranks = list(set(rank_dict.values()))

        new_db = mtdb().set_index()
        for ome in self:
            for rank in ranks:
                try:
                    if self[ome]["taxonomy"][rank].lower() in lineages:
                        new_db[ome] = self[ome]
                except KeyError:  # invalid rank key for row
                    pass  # probably should standardize tax jsons period

        return new_db

    def extract_ome(self, omes: Iterable[str], column: str = "ome") -> "mtdb":
        """Extract a list of genome codes (omes) of interest"""
        new_db = mtdb().set_index(column)
        db = self.set_index(column)
        for i in db:
            if i in list(omes):
                new_db[i] = db[i]
        return new_db.set_index()

    def extract_source(self, source: str) -> "mtdb":
        """Extract an MTDB with genomes from a particular source"""
        return mtdb(
            {
                ome: row
                for ome, row in self.items()
                if row["source"].lower() == source.lower()
            },
            index="ome",
        )

    def extract_pub(self) -> "mtdb":
        """Extract only published and usable genomes"""
        new_db = mtdb().set_index()
        for ome, row in self.items():
            if row["published"]:
                new_db[ome] = row
        return new_db


def get_login(ncbi, jgi):

    ncbi_api, jgi_email, jgi_pwd = None, None, None
    print(flush=True)
    if ncbi:
        ncbi_api = getpass.getpass(prompt="NCBI api key (blank if none): ")
    if jgi:
        jgi_email = input("JGI email: ")
        jgi_pwd = getpass.getpass(prompt="JGI password (required): ")
    print(flush=True)

    return ncbi_api, jgi_email, jgi_pwd


# Path to the UNENCRYPTED credential store (see store_login). Kept separate from
# the password-encrypted key (`~/.mycotools/mtdb_key`) so the two never collide.
PLAIN_LOGIN_PATH = "~/.mycotools/mtdb_credentials.json"


def store_login(
    ncbi_api,
    jgi_email,
    jgi_pwd,
    info_path=PLAIN_LOGIN_PATH,
    encrypted_path="~/.mycotools/mtdb_key",
):
    """Store NCBI/JGI credentials WITHOUT a MycotoolsDB password.

    Credentials are written as JSON with owner-only (0600) permissions. Unlike
    `encrypt_pw`, no password is set or required to read them back - trading
    encryption for convenience. Anyone able to read the file can read the JGI
    password in the clear, so the file permissions are its only protection.
    """
    info_path = format_path(info_path)
    Path(info_path).parent.mkdir(parents=True, exist_ok=True)
    data = {
        "ncbi_api": ncbi_api or "",
        "jgi_email": jgi_email or "",
        "jgi_pwd": jgi_pwd or "",
    }
    # create the file with restrictive permissions *before* writing the secret,
    # so the password is never briefly exposed with a broader umask
    fd = os.open(info_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w") as out:
        json.dump(data, out)
    os.chmod(info_path, 0o600)
    logger.warning(
        "Stored credentials UNENCRYPTED at %s (permissions 600). Anyone able to "
        "read this file can read your JGI password.",
        info_path,
    )
    # a password-encrypted key would otherwise take precedence in login_check;
    # remove it so the no-password store is the one that is actually used
    enc = Path(format_path(encrypted_path))
    if enc.is_file():
        enc.unlink()
        logger.info("Removed prior password-encrypted key %s", str(enc))


def read_plain_login(info_path=PLAIN_LOGIN_PATH):
    """Return (ncbi_api, jgi_email, jgi_pwd) from the unencrypted store written
    by `store_login`. Missing fields come back as empty strings."""
    with open(format_path(info_path), "r") as raw:
        data = json.load(raw)
    return (
        data.get("ncbi_api", ""),
        data.get("jgi_email", ""),
        data.get("jgi_pwd", ""),
    )


def encrypt_pw(
    ncbi_api,
    jgi_email,
    jgi_pwd,
    info_path=format_path("~/.mycotools/mtdb_key"),
):
    salt = b"D9\x82\xbfSibW(\xb1q\xeb\xd1\x84\x118"
    from cryptography.fernet import Fernet
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC

    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=salt,
        iterations=480000,
    )
    print(flush=True)
    hash_pwd, hash_check = True, False
    while hash_pwd != hash_check:
        if hash_pwd != True:
            logger.error("passwords do not match")
        hash_pwd = getpass.getpass(prompt="New MycotoolsDB login password: ")
        hash_check = getpass.getpass(prompt="Confirm password: ")

    key = base64.urlsafe_b64encode(kdf.derive(hash_pwd.encode("utf-8")))
    fernet = Fernet(key)
    out_data = ncbi_api + "\t" + jgi_email + "\t" + jgi_pwd
    encrypt_data = fernet.encrypt(out_data.encode("utf-8"))
    with open(format_path(info_path), "wb") as out:
        out.write(encrypt_data)
    # keep credentials in exactly one place: drop any unencrypted store
    plain = Path(format_path(PLAIN_LOGIN_PATH))
    if plain.is_file():
        plain.unlink()
        logger.info("Removed unencrypted credential store %s", str(plain))


def login_check(info_path="~/.mycotools/mtdb_key", ncbi=True, jgi=True, encrypt=False):
    salt = b"D9\x82\xbfSibW(\xb1q\xeb\xd1\x84\x118"
    # Credential source precedence:
    #   1. password-encrypted key   (encrypt_pw)      - prompts for a password
    #   2. unencrypted store        (store_login)     - no password required
    #   3. interactive prompt       (get_login)        - not persisted
    if Path(format_path(info_path)).is_file():
        from cryptography.fernet import Fernet
        from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
        from cryptography.hazmat.primitives import hashes

        kdf = PBKDF2HMAC(
            algorithm=hashes.SHA256(),
            length=32,
            salt=salt,
            iterations=480000,
        )
        print(flush=True)
        if sys.stdin.isatty():
            hash_pwd = getpass.getpass(prompt="MTDB login password (stdin allowed): ")
        else:
            hash_pwd = sys.stdin.readline().rstrip()
        key = base64.urlsafe_b64encode(kdf.derive(hash_pwd.encode("utf-8")))
        fernet = Fernet(key)
        #        with open(format_path(info_path) + '/.key', 'rb') as raw_key:
        #           fernet = Fernet(raw_key)
        with open(format_path(info_path), "rb") as raw_file:
            data = raw_file.read()
        decrypted = fernet.decrypt(data)
        data = decrypted.decode("UTF-8").split("\t")
        # legacy key files stored a leading NCBI email; drop it if present
        if len(data) == 4:
            data = data[1:]
        if len(data) != 3:
            logger.error("BAD PASSWORD FILE. Delete ~/.mycotools/mtdb_key to reset.")
            sys.exit(8)
        ncbi_api = data[0].rstrip()
        jgi_email = data[1].rstrip()
        jgi_pwd = data[2].rstrip()
        return ncbi_api, jgi_email, jgi_pwd
    elif Path(format_path(PLAIN_LOGIN_PATH)).is_file():
        # unencrypted store written by store_login - no password required
        return read_plain_login(PLAIN_LOGIN_PATH)
    else:
        ncbi_api, jgi_email, jgi_pwd = get_login(ncbi, jgi)
        # CURRENTLY THE REST DOESNT WORK, SO SKIP FOR NOW
        return ncbi_api, jgi_email, jgi_pwd


# opens a `log` file path to read, searches for the `ome` code followed by a whitespace character, and edits the line with `edit`
def log_editor(log, ome, edit):

    try:
        with open(log, "r") as raw:
            data = raw.read()
    except FileNotFoundError:
        data = ""

    if re.search(ome + r"\s", data):
        new_data = re.sub(ome + r"\s.*", edit, data)
    else:
        data = data.rstrip()
        data += "\n"
        new_data = data + edit

    with open(log, "w") as towrite:
        towrite.write(new_data)


def read_log(log, columns="", sep="\t"):

    log_dict = {}
    with open(log, "r") as raw:
        data = raw.read()
    dataLines = [x.split(sep) for x in data.split("\n")]

    if type(columns) is str:
        if dataLines[0][0].startswith("#"):
            dataLines[0][0] = re.sub(r"^\#+", "", dataLines[0][0])
        columns = [head for head in dataLines[0]]
        del dataLines[0]
    elif type(columns) is int:
        x = 0
        columns = [x + 1 for y in dataLines[0]]

    for line in dataLines:
        log_dict[line[0]] = {columns[x]: line[x] for x in range(1, len(line))}

    return log_dict


def primary_db(path="$MYCODB", verbose=True):
    """Acquire the path of the primary database by searching $MYCODB for a file
    with a basename that starts with a date string %Y%m%d."""

    path = path.replace("$", "")
    try:
        full_path = os.environ[path]
    except KeyError:  # $MYCODB not initialized
        return None

    files = collect_files(full_path, "mtdb")
    basenames = [Path(x).name for x in files]
    dates = [x.replace(".mtdb", "") for x in basenames if re.search(r"^\d+\.mtdb$", x)]
    primary = "19991231"  # arbitrary primary for sorting
    for date in dates:
        dtDate = datetime.datetime.strptime(date, "%Y%m%d")
        dtMast = datetime.datetime.strptime(primary, "%Y%m%d")
        if dtDate > dtMast:
            primary = date
    if primary == "19991231":  # if it is the arbitrary start
        if verbose:
            logger.warning("Primary MTDB not found in " + full_path)
        return None
    primary_path = format_path("$" + path + "/" + primary + ".mtdb")

    return primary_path


# imports database, converts into df
# returns database dataframe
def db2df(data, stdin=False):
    """Deprecated legacy Pandas implementation of MTDB import"""
    import pandas as pd

    columns = mtdb.columns
    if isinstance(data, mtdb):
        db_df = pd.DataFrame(data.reset_index())
    elif not stdin:
        data = format_path(data)
        db_df = pd.read_csv(data, sep="\t")
        if "ome" not in set(db_df.columns) and "assembly_acc" not in set(db_df.columns):
            db_df = pd.read_csv(data, sep="\t", header=None)
    else:
        db_df = pd.read_csv(StringIO(data), sep="\t")
        if "ome" not in set(db_df.columns) and "assembly_acc" not in set(db_df.columns):
            db_df = pd.read_csv(StringIO(data), sep="\t", header=None)

    db_df = db_df.fillna("")

    db_df.columns = columns
    for i, row in db_df.iterrows():
        db_df.at[i, "taxonomy"] = read_tax(row["taxonomy"])
        db_df.at[i, "taxonomy"]["genus"] = row["genus"]
        db_df.at[i, "taxonomy"]["species"] = row["genus"] + " " + row["species"]
        # if malformatted due to decreased entries in some lines, this will raise an IndexError
        if (
            row["fna"] or row["fna"] == row["ome"] + ".fna"
        ):  # abbreviated line w/o file coordinates
            db_df.at[i, "fna"] = os.environ["MYCOFNA"] + row["ome"] + ".fna"
            db_df.at[i, "faa"] = os.environ["MYCOFAA"] + row["ome"] + ".faa"
            db_df.at[i, "gff3"] = os.environ["MYCOGFF3"] + row["ome"] + ".gff3"
        else:  # has file coordinates
            db_df.at[i, "fna"] = format_path(row["fna"])
            db_df.at[i, "faa"] = format_path(row["faa"])
            db_df.at[i, "gff3"] = format_path(row["gff3"])

    return db_df


def df2std(df):
    """
    Standardized organization of database columns. DEPRECATED pandas MTDB
    """
    # temporary check, will remove when all scripts account for new biosample column
    trans_df = df[mtdb.columns]
    return trans_df


# imports dataframe and organizes it in accord with the `std df` format
# exports as database file into `db_path`. If rescue is set to 1, save in home folder if output doesn't exist
# if rescue is set to 0, do not output database if output dir does not exit
def df2db(df, db_path, header=False, overwrite=False, std_col=True, rescue=True):
    """Deprecated output pandas MTDB implementation to file"""
    import pandas as pd

    df = df.set_index("ome")
    df = df.sort_index()
    df = df.reset_index()

    if db_path != sys.stdout:
        db_path = format_path(db_path)
    elif overwrite:
        number = 0
        while Path(db_path).exists():
            number += 1
            db_path = str(Path(db_path)) + "_" + str(number)

    if std_col:
        df = df2std(df)

    while True:
        try:
            df.to_csv(db_path, sep="\t", index=None, header=header)
            break
        except FileNotFoundError:
            if rescue:
                logger.warning(
                    "Output directory does not exist. Attempting to save in home folder."
                )
                db_path = "~/" + Path(str(Path(db_path))).name
                df.to_csv(db_path, sep="\t", index=None)
                raise FileNotFoundError
                break
            else:
                logger.error("Output directory does not exist. Rescue not enabled.")
                raise FileNotFoundError
                break


def hit2taxonomy(
    taxid, rank="kingdom", lineage="fungi", skip=False, api=None
):
    """
    Takes a searchTerm string, queries NCBI via Entrez, obtains TaxIDs,
    if there are multiple hits return an error an query the TaxID for
    taxonomy to find a fungal hit. If there is a fungal hit, create a
    tax_dict and return it.
    tax_dict = { taxonomicRank: rankName }
    """

    tax_dict, count, sleep = {}, 0, False
    while True:
        try:
            tax_handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
            break
        except urllib.error.HTTPError as goon:
            count += 1
            if count == 5 and skip:
                logger.error("5 failed HTTP queries. Is NCBI down?")
                sys.exit(100)
            if 500 <= goon.code <= 599:
                time.sleep(1)
            sleep = True

    count = 0
    while True:
        try:
            records = Entrez.read(tax_handle, validate=False)
            break
        except Entrez.Parser.NotXMLError:
            time.sleep(1)
            count += 1
            if count == 5:
                logger.error("5 failed taxonomy queries. " + str(taxid))
                logger.debug("%s", tax_handle)
                for line in tax_handle:
                    logger.debug(line)
                if skip:
                    sys.exit(1)
                records = False
                break
            tax_handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
            if "latin-1" in str(tax_handle):
                count = 0

                logger.error("latin-1 encoding")
                for line in tax_handle:
                    logger.debug(line)
                time.sleep(30)
                Entrez.api_key = api
            sleep = True

    if records:
        if "LineageEx" in records[0]:
            lineages = records[0]["LineageEx"]
            for tax in lineages:
                if tax["Rank"].lower() == rank.lower():
                    tax_dict = {
                        str(x["Rank"]): str(x["ScientificName"]) for x in lineages
                    }
                    continue

    return tax_dict, sleep


def prepare_tax_dicts(df, tax_dicts=None):
    """Identify the genera that do not have higher taxonomy ascribed to them.

    Backwards-compatible dispatcher: MTDBs (and MTDB-shaped dicts) use
    ``mtdb.prepare_tax_dicts``; a deprecated pandas DataFrame uses the path
    below."""
    if tax_dicts is None:
        tax_dicts = {}
    # is the df of mtdb class or have a taxonomy column as a list?
    if isinstance(df, mtdb):
        return df.prepare_tax_dicts(tax_dicts)
    if isinstance(df["taxonomy"], list):
        return mtdb(df).prepare_tax_dicts(tax_dicts)
    # otherwise it is a pandas dataframe
    need_tax = set()
    df["taxonomy"] = df["taxonomy"].fillna({})
    for k, v in df.iterrows():
        if v["genus"] in tax_dicts:
            continue
        tax_json = read_tax(v["taxonomy"])
        if any(
            name
            for rank, name in tax_json.items()
            if rank not in {"genus", "species", "strain"}
        ):
            tax_dicts[v["genus"]] = tax_json
        else:
            need_tax.add(v["genus"])
    need_tax = set(need_tax).difference(set(tax_dicts.keys()))
    return need_tax, tax_dicts


def query_ncbi4taxonomy(genus, api_key, king, rank, count=0):
    """gather taxonomy from information present in the genus column
    and query NCBI using Entrez"""

    # acquire the TaxIDs related to the genus
    ids = []
    for attempt in range(5):
        try:
            search_term = str(genus) + " [GENUS]"
            handle = Entrez.esearch(db="Taxonomy", term=search_term)
            ids = Entrez.read(handle)["IdList"]
            break
        except:
            time.sleep(3)
            count = 0

    count += 1
    if not ids:
        logger.warning(f"{genus} TaxID acquisition failed")
        return None, count

    # for each taxID acquired, fetch the actual taxonomy information
    taxid = 0
    for tax in ids:
        if not api_key and count >= 2:
            time.sleep(1)
            count = 0
        elif api_key and count >= 9:
            time.sleep(1)
            count = 0

    for attempt in range(5):
        try:
            handle = Entrez.efetch(db="Taxonomy", id=tax, remode="xml")
            records = Entrez.read(handle)
            if records:
                lineages = records[0]["LineageEx"]
            else:
                lineages = []
            break
        except:
            time.sleep(3)
            lineages = []
            count = 0
    #            handle = Entrez.efetch(db="Taxonomy", id=tax, remode = "xml")
    #           records = Entrez.read(handle)
    #          lineages = records[0]['LineageEx']
    #         break

    count += 1
    # for each lineage, if it is a part of the kingdom, `king`,
    # use that TaxID if there are multiple TaxIDs, use the first one found
    if king:
        for lineage in lineages:
            if lineage["Rank"].lower() == rank.lower():
                if lineage["ScientificName"].lower() == king.lower():
                    taxid = tax
                    if len(ids) > 1:
                        logger.warning("Multiple Tax IDs: " + str(ids))
                    break
    else:
        for lineage in lineages:
            taxid = tax

    if taxid == 0:
        logger.warning(f"{genus} not recovered in {king}")
        return None, count

    # for each taxonomic classification, add it to the taxonomy dictionary string
    # append each taxonomy dictionary string to a list of dicts
    genus_tax = {}
    for tax in lineages:
        genus_tax[tax["Rank"].lower()] = tax["ScientificName"]

    return genus_tax, count


def gather_taxonomy(
    df, api_key=None, king="fungi", rank="kingdom", tax_dicts={}, output_path=None
):
    """Gather taxonomy for an MTDB or deprecated pd MTDB by querying NCBI's
    taxonomy hierarchy"""
    need_tax, tax_dicts = prepare_tax_dicts(df, tax_dicts)

    count = 0
    for genus in tqdm(sorted(need_tax), total=len(need_tax)):
        # if there is no api key, sleep for a second after 3 queries
        # if there is an api key, sleep for a second after 10 queries
        if not api_key and count >= 2:
            time.sleep(1)
            count = 0
        elif api_key and count >= 9:
            time.sleep(1)
            count = 0

        genus_tax, count = query_ncbi4taxonomy(genus, api_key, king, rank, count)
        if genus_tax:
            tax_dicts[genus] = genus_tax

    if output_path:
        with open(output_path + ".tmp", "w") as out:
            for genus, tax_dict in tax_dicts.items():
                out.write(f"{genus}\t{json.dumps(tax_dict)}\n")
        Path(f"{output_path}.tmp").rename(output_path)
    return tax_dicts


def gather_taxonomy_dataset(
    df,
    api_key=None,
    king="fungi",
    rank="kingdom",
    tax_dicts={},
    output_path=None,
    verbose=False,
):
    """Gather taxonomy for an MTDB or deprecated pd MTDB by querying NCBI's
    taxonomy hierarchy"""
    if not output_path:
        tax_dicts = gather_taxonomy_fallback(
            df, api_key=api_key, king=key, rank=rank, tax_dicts=tax_dicts
        )
        return tax_dicts

    need_tax, tax_dicts = prepare_tax_dicts(df, tax_dicts)

    tax_accs_file = output_path + "query_genera.txt"
    with open(tax_accs_file, "w") as out:
        out.write("\n".join(sorted(need_tax)))

    cwd = str(Path.cwd())
    os.chdir(output_path)
    cmd_scaf = [
        "datasets",
        "download",
        "taxonomy",
        "taxon",
        "--inputfile",
        tax_accs_file,
    ]
    if api_key:
        cmd_scaf.extend(["--api-key", api_key])

    #    if verbose:
    #       v = None
    #  else:
    #     v = subprocess.DEVNULL
    v = None
    count = 0
    dataset_path = output_path + "ncbi_dataset.zip"
    while count < 3:
        cmd_call = subprocess.call(cmd_scaf, stdout=v, stderr=v)
        try:
            with zipfile.ZipFile(dataset_path, "r") as zip_ref:
                zip_ref.extractall(zip_ref)
        except zipfile.BadZipFile:
            count += 1
            if count == 3:
                logger.error("taxonomy acquisition failed")
                sys.exit(130)
            continue
        break

    Path(dataset_path).unlink()
    unzip_path = output_path + "ncbi_dataset/"

    tax_dicts = parse_dataset_taxonomy(
        f"{unzip_path}data/taxonomy_report.jsonl", tax_dicts, king, rank
    )

    os.chdir(cwd)
    return tax_dicts


def parse_dataset_taxonomy(tax_json, tax_dicts, tax_head, rank_head):
    """Parse the dataset taxonomy list"""
    tax_strs = {
        "superkingdom",
        "kingdom",
        "phylum",
        "subphylum",
        "class",
        "order",
        "family",
        "subfamily",
    }

    with open(tax_json, "r") as raw:
        for line in raw:
            tax_line = json.loads(line.rstrip())
            tax_dict = {}
            if rank_head.lower() in tax_line["classification"]:
                if (
                    tax_line["classification"][rank_head.lower()].lower()
                    != tax_head.lower()
                ):
                    continue
                for rank, data in tax_line["classification"]:
                    if rank.lower() in tax_strs:
                        tax = data["name"]
                        tax_dict[rank] = tax
                genus = tax_dict["genus"]["name"]
                tax_dicts[genus] = tax_dict

    return tax_dicts


# read_tax is defined as mtdb.read_tax; expose it at module scope for the
# historical ``from mycotools.lib.dbtools import read_tax`` import.
read_tax = mtdb.read_tax


# assimilate taxonomy dictionary strings and append the resulting taxonomy string dicts to an inputted database
# forbid a list of taxonomic classifications you are not interested in and return a new database
def assimilate_tax(db, tax_dicts, ome_index="ome", forbid=None):
    """Backwards-compatible dispatcher for ``mtdb.assimilate_tax``; retains the
    deprecated pandas DataFrame path."""
    if forbid is None:
        forbid = mtdb._forbidden_tax_ranks
    if isinstance(db, mtdb):
        return db.assimilate_tax(tax_dicts, forbid=forbid)
    tax_dicts = mtdb._reconcile_tax_dicts(set(db["genus"]), tax_dicts, forbid)
    for i, row in db.iterrows():
        db.at[i, "taxonomy"] = tax_dicts[row["genus"]]
    return db, tax_dicts


def parse_user_config(mtdb_config_file=format_path("~/.mycotools/config.json")):
    config_dir = format_path("~/.mycotools/")
    if not Path(config_dir).is_dir():
        Path(config_dir).mkdir()
        config_dir += "/"

    if Path(mtdb_config_file).is_file():
        config = read_json(mtdb_config_file)
    else:
        config = {"log": {}}
    return config


def mtdb_disconnect(mtdb_config_file=format_path("~/.mycotools/config.json")):
    user_config = parse_user_config()
    user_config["active"] = False
    write_json(user_config, mtdb_config_file)


def mtdb_connect(
    user_config, mycodb_loc, mtdb_config_file=format_path("~/.mycotools/config.json")
):
    user_config["active"] = mycodb_loc
    write_json(user_config, mtdb_config_file)
    for var, env in user_config[user_config["active"]].items():
        os.environ[var] = env


def mtdb_initialize(
    mycodb_loc, mtdb_config_file=format_path("~/.mycotools/config.json"), init=False
):

    user_config = parse_user_config()
    # to maintain legacy installs
    if "log" not in user_config:
        user_config["log"] = {}

    mtdb_config = read_json(mycodb_loc + "config/mtdb.json")
    dbtype = mtdb_config["branch"]
    logger.info("Establishing " + dbtype + " connection")

    if not Path(mycodb_loc + "mtdb/").is_dir() and not init:
        raise FileNotFoundError("invalid MycotoolsDB path")
    dPath = mycodb_loc + "data/"
    user_config[mycodb_loc] = {
        "MYCODB": mycodb_loc + "mtdb/",
        "MYCOFNA": dPath + "fna/",
        "MYCOFAA": dPath + "faa/",
        "MYCOGFF3": dPath + "gff3/",
    }

    login_time = datetime.datetime.now().strftime("%Y%m%d %H:%M:%S")
    #    login_time = datetime.datetime.now().strftime('%Y%m%d')
    if dbtype in user_config["log"]:
        user_config["log"][dbtype][mycodb_loc] = login_time
    else:
        user_config["log"][dbtype] = {mycodb_loc: login_time}

    mtdb_connect(user_config, mycodb_loc, mtdb_config_file)


interface = format_path("~/.mycotools/config.json")
if Path(interface).is_file():
    envs_info = read_json(interface)
    if envs_info["active"]:
        for var, env in envs_info[envs_info["active"]].items():
            os.environ[var] = env

# if not primary_db():
#    eprint('WARNING: Primary MycotoolsDB not connected; setup using `mtdb u/-i/-p/-f`', flush = True)

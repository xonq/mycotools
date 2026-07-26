#! /usr/bin/env python3

import os
import re
import sys
import glob
import gzip
import json
import logging
import shutil
import tarfile
import argparse
import subprocess
from pathlib import Path
from contextlib import contextmanager
from tqdm import tqdm
from datetime import datetime


logger = logging.getLogger(__name__)


class _LevelFormatter(logging.Formatter):
    """Prefix records with their level name, except INFO, which is emitted
    verbatim so status messages read as they did before the logging
    migration."""

    _PREFIX = {
        logging.DEBUG: "DEBUG: ",
        logging.INFO: "",
        logging.WARNING: "WARNING: ",
        logging.ERROR: "ERROR: ",
        logging.CRITICAL: "CRITICAL: ",
    }

    def format(self, record):
        message = super().format(record)
        prefix = self._PREFIX.get(record.levelno, "")
        # Do not double-prefix messages that already carry their level word.
        if prefix and message.lstrip().upper().startswith(prefix.strip().rstrip(":")):
            prefix = ""
        return prefix + message


def setup_logging(verbose=False, level=None):
    """Configure logging for the Mycotools CLI.

    Attaches a single stderr handler to the root logger (idempotent across
    repeated calls) and sets the ``mycotools`` logger to INFO, or DEBUG when
    ``verbose`` is True. The root logger stays at WARNING so third-party
    libraries do not clutter diagnostic output. Call once from a script's
    entry point after parsing arguments."""
    if level is None:
        level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()
    root.setLevel(logging.WARNING)
    if not any(getattr(h, "_mycotools", False) for h in root.handlers):
        handler = logging.StreamHandler(sys.stderr)
#        handler.setFormatter(_LevelFormatter("%(message)s"))
        handler._mycotools = True
        root.addHandler(handler)
    logging.getLogger("mycotools").setLevel(level)
    return logging.getLogger("mycotools")


class KonLog:
    """A print class designed to enable swift string formatting while moving
    between scripts"""

    def __init__(self, base=0, spaces=4, flush=True, tee=None):
        """Initialize a print class with the base index, the number of spaces
        between bases, whether or not to flush after printing, and a
        placeholder for a tee variable that will input the path of a log output
        file to write to"""
        self._base = base
        self._spaces = spaces
        self.flush = flush

    def print(self, print_str, add=0, stderr=False, flush=None, **kwargs):
        """Perform the basic print function"""
        if flush is None:
            flush = self.flush
        else:
            flush = flush
        if stderr:
            print(
                f"{print_str:>{(self._base + add) * self._spaces}}",
                file=sys.stderr,
                **kwargs,
            )
        else:
            print(f"{print_str:>{(self._base + add) * self._spaces}}", **kwargs)

    def up(self):
        self._base += 1

    def down(self):
        if self._base > 0:
            self._base -= 1

    @property
    def base(self):
        return self._base

    @base.setter
    def base(self, value):
        if value >= 0:
            self._x = value
        else:
            raise ValueError("base value for print log must be >= 0")

    @base.deleter
    def base(self):
        self._base = 0


def checksum(path, cmd="sha256", ref=""):
    if not isinstance(path, list):
        hash_res = subprocess.run([cmd + "sum", path], stdout=subprocess.PIPE)
        hash_info = hash_res.stdout.decode("utf-8")
        hash_find = re.search(r"\w+", hash_info)
        hash_ = hash_find[0]
        if ref:
            if hash_ == ref:
                return True
            else:
                return False
        else:
            return hash_
    else:
        path0, path1 = path[0], path[1]
        hash_cmd = subprocess.call([cmd + "sum", path0, path1])
        return bool(not hash_cmd)


def stdin2str():
    data = ""
    for line in sys.stdin:
        data += line.rstrip() + "\n"
    data = data.rstrip()
    return data


def split_input(inp, ignore_white=False):
    """Split a standard input by removing quotations, commas, etc. and
    splitting a list from various bash inputs based on white space"""
    if not inp:
        return []

    if inp.startswith('"') and inp.endswith('"'):
        inp = inp[1:-1]
    elif inp.startswith("'") and inp.endswith("'"):
        inp = inp[1:-1]
    if ", " in inp:
        out_list = inp.split(", ")
    elif "; " in inp:
        out_list = inp.split("; ")
    elif "," in inp:
        out_list = inp.split(",")
    elif not ignore_white:
        out_list = inp.split()
    return out_list


def namespace_to_dict(namespace):
    return {
        k: namespace_to_dict(v) if isinstance(v, argparse.Namespace) else v
        for k, v in vars(namespace).items()
    }


def parse_run_log(log_path, args_dict, fail=set()):
    """Parse a log of the previous run, or create one if it doesn't exist.
    Raise a KeyError if there are deviations that should fail the existing
    run"""
    if isinstance(args_dict, argparse.Namespace):
        args_dict = namespace_to_dict(args_dict)

    if not Path(log_path).is_file():
        write_json(args_dict, log_path)
        return {}
    else:
        prev_args = read_json(log_path)
        deviations = {}
        for arg, val in args_dict.items():
            if arg in prev_args:
                if val != prev_args[arg]:
                    deviations[arg] = prev_args[arg]
            else:
                deviations[arg] = None
        failure_args = sorted(x for x in deviations if x in fail)
        if failure_args:
            raise KeyError(
                f"arguments discrepant from previous run: " + f"{failure_args}"
            )
        return deviations


def hex2rgb(hexCode):
    return tuple(int(hexCode.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))


def get_colors(size, ignore=[], rgb=False):
    if size < 16:
        colors = [
            "#000000",
            "#004949",
            "#009292",
            "#ff6db6",
            "#ffb6db",
            "#490092",
            "#006ddb",
            "#b66dff",
            "#6db6ff",
            "#b6dbff",
            "#920000",
            "#924900",
            "#db6d00",
            "#24ff24",
            "#ffff6d",
        ]
    elif size < 27:
        colors = [
            "#F0A3FF",
            "#0075DC",
            "#993F00",
            "#4C005C",
            "#191919",
            "#005C31",
            "#2BCE48",
            "#FFCC99",
            "#808080",
            "#94FFB5",
            "#8F7C00",
            "#9DCC00",
            "#C20088",
            "#003380",
            "#FFA405",
            "#FFA8BB",
            "#426600",
            "#FF0010",
            "#5EF1F2",
            "#00998F",
            "#E0FF66",
            "#740AFF",
            "#990000",
            "#FFFF80",
            "#FFFF00",
            "#FF5005",
        ]
    else:
        colors = [
            "#000000",
            "#010067",
            "#d5ff00",
            "#ff0056",
            "#9e008e",
            "#0e4ca1",
            "#ffe502",
            "#005f39",
            "#00ff00",
            "#95003a",
            "#ff937e",
            "#a42400",
            "#001544",
            "#91d0cb",
            "#620e00",
            "#6b6882",
            "#0000ff",
            "#007db5",
            "#6a826c",
            "#00ae7e",
            "#c28c9f",
            "#be9970",
            "#008f9c",
            "#5fad4e",
            "#ff0000",
            "#ff00f6",
            "#ff029d",
            "#683d3b",
            "#ff74a3",
            "#968ae8",
            "#98ff52",
            "#a75740",
            "#01fffe",
            "#ffeee8",
            "#fe8900",
            "#bdc6ff",
            "#01d0ff",
            "#bb8800",
            "#7544b1",
            "#a5ffd2",
            "#ffa6fe",
            "#774d00",
            "#7a4782",
            "#263400",
            "#004754",
            "#43002c",
            "#b500ff",
            "#ffb167",
            "#ffdb66",
            "#90fb92",
            "#7e2dd2",
            "#bdd393",
            "#e56ffe",
            "#deff74",
            "#00ff78",
            "#009bff",
            "#006401",
            "#0076ff",
            "#85a900",
            "#00b917",
            "#788231",
            "#00ffc6",
            "#ff6e41",
            "#e85ebe",
        ]
    if rgb:
        rgbColors = []
        for color in colors:
            rgbColors.append(hex2rgb(color))
        colors = rgbColors

    if ignore:
        if isinstance(ignore, str):
            ignore = [ignore]
        for ig in ignore:
            try:
                colors.pop(colors.index(ig))
            except ValueError:
                pass

    return colors


def tardir(dir_, rm=True):

    if not Path(format_path(dir_)).is_dir():
        return False
    with tarfile.open(format_path(dir_)[:-1] + ".tar.gz", "w:gz") as tar:
        tar.add(dir_, arcname=Path(format_path(dir_)[:-1]).name)
    if rm:
        shutil.rmtree(dir_)


def untardir(dir_, rm=False, to=None):
    if not to:
        to = str(Path(dir_[:-1]).parent)
    tar = tarfile.TarFile.open(dir_)
    tar.extractall(path=to)
    tar.close()
    if rm:
        Path(dir_).unlink()


def checkdir(dir_, unzip=False, to=None, rm=False):
    if Path(dir_).is_dir():
        return True
    elif Path(format_path(dir_) + ".tar.gz").is_file():
        if unzip:
            if dir_.endswith("/"):
                dir_ = dir_[:-1]
            untardir(dir_ + ".tar.gz", rm=rm, to=to)
            return True
    return False


def read_json(config_path, compress=False):

    if compress or config_path.endswith(".gz"):
        with gzip.open(config_path, "rt") as json_raw:
            json_dict = json.load(json_raw)
    else:
        with open(config_path, "r") as json_raw:
            json_dict = json.load(json_raw)

    return json_dict


def write_json(obj, json_path, compress=False, indent=1, **kwargs):
    if compress or json_path.endswith(".gz"):
        with gzip.open(json_path, "wt") as json_out:
            json.dump(obj, json_out, indent=indent, **kwargs)
    else:
        with open(json_path, "w") as json_out:
            json.dump(obj, json_out, indent=indent, **kwargs)


@contextmanager
def atomic_write(path, mode="w", suffix=".tmp", **kwargs):
    """Context manager for crash-safe file writes.

    Yields a handle to a temporary sibling file (`path` + `suffix`) and, only
    upon leaving the block without error, atomically replaces `path` with it.
    Because the temporary file lives beside the target it shares a filesystem,
    so the final replace is atomic and can never leave a half-written `path`
    behind. If the block raises, the temporary file is removed and `path` is
    left untouched.

    Consolidates the `open(path + '.tmp')` ... `shutil.move`/`Path.rename`
    idiom repeated throughout the codebase.
    """
    path = Path(path)
    tmp_path = path.with_name(path.name + suffix)
    try:
        with open(tmp_path, mode, **kwargs) as handle:
            yield handle
        tmp_path.replace(path)
    except BaseException:
        if tmp_path.is_file():
            tmp_path.unlink()
        raise


def gunzip(gzip_file, remove=True, spacer="\t"):
    """gunzips gzip_file and removes if successful"""

    new_file = re.sub(r"\.gz$", "", format_path(gzip_file))
    try:
        with gzip.open(format_path(gzip_file), "rt") as f_in:
            with open(new_file, "w") as f_out:
                for line in f_in:
                    f_out.write(line)
        if remove:
            Path(gzip_file).unlink()
        return new_file
    except:
        if Path(new_file).is_file():
            if Path(gzip_file).is_file():
                Path(new_file).unlink()
        raise IOError("gunzip " + str(gzip_file) + " failed")


def fmt_float(val, sig_dig=None):

    val_str = str(val)
    try:
        e = val_str.index("e")
        dig = val_str[:e].replace(".", "")
        num = val_str[e + 1 :]
        if num.startswith("-"):
            val_op = "0."
            for i in range(abs(int(num)) - 1):
                val_op += "0"
            val_str = val_op + dig.replace(".", "")
        else:
            num = abs(int(num))
            if len(dig) > num:
                val_op = list(dig)
                val_op.insert(".", num + 1)
                val_str = "".join(val_op)
            else:
                zeroes = num + 1 - len(dig)
                val_op = dig
                for i in range(zeroes):
                    val_op += "0"
                val_str = val_op
    except ValueError:
        val_str = str(val)

    if sig_dig:
        if len(val_str.replace(".", "")) > sig_dig:
            try:
                per_i = val_str.index(".")
            except ValueError:
                per_i = None
            val_list = list(val_str)
            if per_i:
                if sig_dig > per_i:
                    i = val_list[sig_dig]
                    post = int(val_list[sig_dig + 1])
                    if post >= 5:
                        val_list.insert(sig_dig, int(i) + 1)
                        val_list.pop(sig_dig + 1)
                    val_str = "".join(str(v) for v in val_list[: sig_dig + 1])
                else:
                    i = val_list[sig_dig - 1]
                    try:
                        post = val_list[sig_dig]
                    except IndexError:
                        post = 0
                    if post == ".":
                        post = val_list[sig_dig + 1]
                    elif i == ".":
                        i = val_list[sig_dig]
                        try:
                            post = val_list[sig_dig + 1]
                        except IndexError:
                            post = 0
                    if int(post) >= 5:
                        val_list.insert(sig_dig, int(i) + 1)
                        val_list.pop(sig_dig + 1)
                    val_list = val_list[:sig_dig]
                    while len(val_str) > len(val_list):
                        val_list.append("0")
                    val_str = "".join(str(v) for v in val_list[: sig_dig + 3])
            else:
                i = val_list[sig_dig + 1]
                post = int(val_list[sig_dig + 2])
                if post == ".":
                    post = int(val_list[sig_dig + 3])
                elif i == ".":
                    i = val_list[sig_dig + 2]
                    try:
                        post = val_list[sig_dig + 3]
                    except IndexError:
                        post = 0
                if int(post) >= 5:
                    val_list.insert(sig_dig + 1, int(i) + 1)
                    val_list.pop(sig_dig + 2)
                val_list = val_list[: sig_dig + 2]
                while len(val_str) > len(val_list):
                    val_list.append("0")
                val_str = "".join(str(v) for v in val_list[: sig_dig + 1])

    return val_str


def find_execs(deps, exit=set(), verbose=True):
    """
    Inputs list of dependencies, `dep`, to check path.
    If dependency is in exit and dependency is not in path,
    then exit.
    """

    logger.debug("Dependency check:")
    checks, failed = [], []
    if isinstance(deps, str):
        deps = [deps]
    for dep in sorted(deps):
        check = shutil.which(dep)
        logger.debug("{:<15}".format(dep + ":") + str(check))
        if not check and dep in exit:
            failed.append(dep)
        else:
            checks.append(check)

    if failed:
        logger.error("missing dependencies: " + ", ".join(failed))
        sys.exit(300)

    return checks


def find_envs(envs, exit=set(), verbose=True):
    """
    Inputs list of paths, `envs`, to check path.
    If env is not in path and it is in exit, exit.
    """

    logger.debug("Environment check:")
    if type(envs) is str:
        envs = [envs]
    for env in envs:
        try:
            logger.debug("{:<15}".format(env + ":") + str(os.environ[env]))
        except KeyError:
            logger.debug("{:<15}".format(env + ":") + "None")
            if env in exit:
                logger.error(env + " not in PATH")
                sys.exit(301)


def format_path(path, force_dir=False):
    """Convert a path to an absolute path with explicit directories.

    Expands ``~`` and ``$VAR`` (raising KeyError on an undefined variable). By
    convention an existing directory is returned with a trailing ``/`` while
    files and non-existent paths have none; ``force_dir`` forces a trailing
    ``/`` for directories that do not exist yet. Symlinks are not resolved."""

    if path:
        path = str(Path(str(path)).expanduser())
        for env in re.findall(r"\$[^/]+", path):
            path = path.replace(env, os.environ[env.replace("$", "")])
        path = path.replace("//", "/")

        if force_dir:
            if not path.endswith("/"):
                path += "/"
        else:
            if path.endswith("/"):
                if not Path(path).is_dir():
                    path = path[:-1]
            else:
                if Path(path).is_dir():
                    path += "/"
            if not path.startswith("/"):
                path = str(Path.cwd()) + "/" + path

        path = path.replace("/./", "/")
        while "/../" in path:
            path = re.sub(r"[^/]+/\.\./(.*)", r"\1", path)

    return path


# need to change recursive to false
def collect_files(directory="./", filetype="*", recursive=False):
    """
    Inputs: directory path, file extension (no "."), recursivity bool
    Outputs: list of files with `filetype`
    If the filetype is a list, split it, else make a list of the one entry.
    Parse the environment variable if applicable. Then, obtain a clean, full
    version of the input directory. Glob to obtain the filelist for each
    filetype based on whether or not it is recursive.
    """

    if type(filetype) == list:
        filetypes = filetype.split()
    else:
        filetypes = [filetype]

    directory = Path(format_path(directory))
    filelist = []
    for filetype in filetypes:
        pattern = "*." + filetype
        matches = directory.rglob(pattern) if recursive else directory.glob(pattern)
        filelist.extend(str(match) for match in matches)

    return filelist


def collect_dirs(input_glob, recursive=False):
    """
    Inputs: directory path, recursive search boolean
    Outputs: list of folders
    Get folders via os, recursively extend output list if True.
    """

    in_dirs = glob.glob(input_glob, recursive=recursive)
    return [dir_ for dir_ in in_dirs if Path(dir_).is_dir()]


def dict_split(Dict, factor):
    """
    Inputs: a dictionary `Dict`, and an integer `factor` to split by
    Outputs: a list of split dictionaries `list_dict`
    Split dictionary keys into a numpy array via the factor. For each
    key, append a new dictionary to the output list and populate it
    with the list of keys within that factor. The last entry will be
    whatever amount is leftover.
    """

    list_dict = []
    keys_list = np.array_split(list(Dict.keys()), factor)

    try:
        for keys in keys_list:
            list_dict.append({})
            for key in keys:
                list_dict[-1][key] = Dict[key]
    except TypeError:
        list_dict = []
        for keys in keys_list:
            list_dict.append({})
            for key in keys:
                list_dict[-1][tuple(key)] = Dict[tuple(key)]

    return list_dict


def sys_start(args, usage, min_len, dirs=[], files=[]):
    """args is a sys.argv typically, or a list of arguments.
    usage is the usage statement without formatting.
    min_len is the minimum number of arguments that should exist.
    dirs is a list of items that should be directories
    files is a list of items that should be files"""

    if "-h" in args or "--help" in args:
        print("\n" + usage + "\n", flush=True)
        sys.exit(0)
    elif len(args) < min_len:
        print("\n" + usage + "\n", flush=True)
        sys.exit(1)
    elif not all(Path(format_path(x)).is_file() for x in files):
        print("\n" + usage, flush=True)
        logger.error("input file(s) do not exist")
        sys.exit(3)
    elif not all(Path(format_path(x)).is_file() for x in dirs):
        print("\n" + usage, flush=True)
        logger.error("input directory does not exist")
        sys.exit(4)

    return args


def inject_args(args, injection_calls):
    manual_cmds = []
    for in_call in injection_calls:
        if in_call in args:
            man_index = args.index(in_call)
            for char in args[man_index + 1]:
                if char == '"' or char == "'":
                    quote_char = char

            manual_cmd, start = [], False
            for i, v in args[man_index + 1 :]:
                if quote_char in v:
                    if not start:
                        start = True
                        manual_cmd.append(v[v.find(quote_char) :])
                    else:
                        manual_cmd.append(v[: v.find(quote_char)])
                        break
                else:
                    manual_cmd.append(v)
            for index in reversed(range(man_index, i + 1)):
                del args[index]
        else:
            manual_cmd = []
        manual_cmds.append(manual_cmd)

    return args, manual_cmds


def intro(script_name, args_dict, credit="", log=False, stdout=True):
    """
    Inputs: script_name string, args_dict dictionary of arguments,
    credit string bool / path for output log
    Outputs: prints an introduction, returns start_time in YYYYmmdd format
    Creates a string to populate and format for the introduction using
    keys as the left-most descriptor and arguments (values) as the right
    most. Optionally outputs a log according to `log` path.
    """

    start_time = datetime.now()

    out_str = (
        "\n" + script_name + "\n" + credit + "\nExecution began: " + str(start_time)
    )

    for arg in args_dict:
        out_str += "\n" + "{:<30}".format(arg.upper() + ":") + str(args_dict[arg])

    logger.info(out_str)

    return start_time


def outro(start_time, log=False, stdout=True):
    """
    Inputs: start time string formatted YYYYmmdd, log path
    Outputs: prints execution time and exits with 0 status
    """

    end_time = datetime.now()
    duration = end_time - start_time
    dur_min = duration.seconds / 60
    out_str = (
        "\nExecution finished: "
        + str(end_time)
        + "\t"
        + "\n\t{:.2}".format(dur_min)
        + " minutes\n"
    )

    logger.info(out_str)

    sys.exit(0)


def prep_output(output, mkdir=True, require_newdir=False, cd=False):
    """
    Inputs: output path, bool `mkdir`, bool `require_newdir`, bool cd.
    Outputs: creates directory, returns None if bool inputs fail, returns
    formatted path
    """

    output = format_path(output, force_dir=True)
    if Path(output).is_dir():
        if require_newdir:
            logger.error("directory exists.")
            return None
    elif Path(output).exists():
        output = str(Path(output).parent)
    else:
        if not mkdir:
            logger.error("directory does not exist.")
            return None
        Path(output).mkdir()

    if cd:
        os.chdir(output)

    return output


def mk_output(base_dir, program, reuse=True, suffix=datetime.now().strftime("%Y%m%d")):
    if not base_dir:
        base_dir = str(Path.cwd()) + "/"
    if not Path(format_path(base_dir)).is_dir():
        raise FileNotFoundError(base_dir + " does not exist")
    if suffix:
        out_dir = format_path(base_dir) + program + "_" + suffix
    else:
        out_dir = format_path(base_dir) + program

    if not reuse:
        count, count_dir = 1, out_dir
        while Path(count_dir).is_dir():
            count_dir += "_" + str(count)
            count += 1
        Path(count_dir).mkdir()
        return count_dir + "/"
    else:
        if not Path(out_dir).is_dir():
            Path(out_dir).mkdir()
        return out_dir + "/"


def check_dep(dep_list=[], var_list=[], exempt=set()):
    """Checks all dependencies in path from list, optional exemption set"""

    failedVars = []
    check = [shutil.which(dep) for dep in dep_list if dep not in exempt]
    for var in var_list:
        if var not in exempt:
            try:
                os.environ[var]
                check.append(True)
            except KeyError:
                check.append(False)
                failedVars.append(var)
    if not all(check):
        logger.error("Dependencies not met:")
        for dep in dep_list:
            if not shutil.which(dep):
                logger.error(dep + " not in PATH")
        for failed in failedVars:
            logger.error(failed + " variable not set")
        sys.exit(135)


def file2list(file_, types="", sep=None, col=None, compress=False):
    """
    Inputs: `file_` path, separator string, column integer index
    Outputs: reads file and returns a list given arguments
    If the `sep` is not valid, return None, else read the file and
    split each line, if `col`, split each by the delimiter, and
    acquire the column index at each line. If not `col`, then
    simply split by the separator and grab the only entry.
    """

    if not compress:
        with open(file_, "r") as raw_data:
            data = raw_data.read() + "\n"
    else:
        with gzip.open(file_, "rt") as raw_data:
            data = raw_data.read() + "\n"

    if col:
        check = data.split("\n")
        if sep:
            check_col = check[0].split(sep).index(col)
            data1_list = [
                x[check_col].rstrip()
                for x in [y.split(sep) for y in check[1:]]
                if len(x) >= check_col + 1
            ]
        else:
            check_col = check[0].split().index(col)
            data1_list = [
                x[check_col].rstrip()
                for x in [y.split() for y in check[1:]]
                if len(x) >= check_col + 1
            ]
    elif types:
        if sep:
            data_list = data.split(sep=sep)
        else:
            data_list = data.split()
        if types == "int":
            try:
                data1_list = [int(x.rstrip()) for x in data_list if x]
            except ValueError:
                data1_list = [int(float(x.rstrip())) for x in data_list if x]
        elif types == "float":
            data_list = [float(x.rstrip()) for x in data_list if x]
    else:
        if sep:
            data_list = data.split(sep=sep)
        else:
            data_list = data.split()
        data1_list = [x.rstrip() for x in data_list if x]

    return data1_list


class Subqueue:
    """Run a subqueue of args; can run commands consecutively to avoid shell
    limitations while maintaining shell injection security"""

    def __init__(self, args, shell=False, verbose=False, injectable=False):
        self.injectable = injectable
        if verbose == 1 or verbose == True:
            self.stdout = None
            self.stderr = None
        elif verbose == 2:
            self.stdout = subprocess.DEVNULL
            self.stderr = None
        else:
            self.stdout = subprocess.DEVNULL
            self.stderr = subprocess.DEVNULL
        if isinstance(args, list) or isinstance(args, tuple) or isinstance(args, set):
            if isinstance(list(args)[0], str):
                self.args = (tuple(args),)
            elif (
                isinstance(list(args)[0], list)
                or isinstance(list(args)[0], tuple)
                or isinstance(list(args)[0], set)
            ):
                self.args = tuple([tuple(x) for x in args])
            else:
                raise TypeError("invalid argument subtype: " + str(type(list(args[0]))))

        elif isinstance(args, str):
            self.args = ((args,),)
        else:
            raise TypeError("invalid arguments type: " + str(type(args)))
        self.complete = False
        self.shell = shell
        self.status = [-1 for x in args]
        self.exit = -1

    def open(self):
        for i, arg in enumerate(self.args):
            if arg[-1] == "&&" and self.injectable:
                failstop = True
                arg = arg[:-1]
            else:
                failstop = False
            self.status[i] = subprocess.call(
                arg, shell=self.shell, stdout=self.stdout, stderr=self.stderr
            )
            if failstop and self.status[i]:
                break
        self.exit = any(x for x in self.status)
        self.complete = True


def run_subqueue(i, args, shell=False, verbose=False, injectable=False):
    #   if args:
    s = Subqueue(args, shell=shell, verbose=verbose, injectable=injectable)
    s.open()
    return i, s.exit


#    return [], -420


def multisub(
    args_lists, processes=1, shell=False, verbose=False, injectable=False, status=True
):
    import multiprocessing as mp

    """
    Inputs: list of arguments, integer of processes, subprocess `shell` bool
    Outputs: launches and monitors subprocesses and returns list of exit 
    information
    If processes are too few, opt to default. While there are entries in the
    `args_lists`: 
        while there are less entries in `running` than max processes
    and arguments remaining in the `args_lists`: Popen the first argument in 
    the list, append the subprocess information and command to `running` and
    delete the argument from `args_lists`. 
        If there are still args in the args_lists, for each index in running
        grab the handle, poll if it is done, and if it is done output the 
        return code and stdin information to outputs and delete the index
        from `running`.
        IF there are not arguments in args_lists, wait until all processes
        are complete before exiting the function.
    """
    if not status:
        with mp.Pool(processes=processes) as pool:
            exit_pre = pool.starmap(
                run_subqueue,
                ((i, x, shell, verbose, injectable) for i, x in enumerate(args_lists)),
            )
    else:
        with mp.Pool(processes=processes) as pool:
            exit_pre = pool.starmap(
                run_subqueue,
                tqdm(
                    [
                        (i, x, shell, verbose, injectable)
                        for i, x in enumerate(args_lists)
                    ],
                    total=len(args_lists),
                ),
            )
    return [exit for i, exit in sorted(exit_pre, key=lambda x: x[0])]

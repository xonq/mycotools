#! /usr/bin/env python3

# NEED to convert gff list to appropriate types

import re
from collections import defaultdict
from itertools import chain

aa_weights = {
    "A": 89.1,
    "R": 174.2,
    "N": 132.1,
    "D": 133.1,
    "C": 121.2,
    "E": 147.1,
    "Q": 146.2,
    "G": 75.1,
    "H": 155.2,
    "I": 131.2,
    "L": 131.2,
    "K": 146.2,
    "M": 149.2,
    "F": 165.2,
    "P": 115.1,
    "S": 105.1,
    "T": 119.1,
    "W": 204.2,
    "Y": 181.2,
    "V": 117.1,
}
aa_weights["X"] = sum(aa_weights.values()) / len(aa_weights)

# codon translation table
translation_table = {
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "AAT": "N",
    "AAC": "N",
    "GAT": "D",
    "GAC": "D",
    "TGT": "C",
    "TGC": "C",
    "CAA": "Q",
    "CAG": "Q",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "CAT": "H",
    "CAC": "H",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "AAA": "K",
    "AAG": "K",
    "ATG": "M",
    "TTT": "F",
    "TTC": "F",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "TGG": "W",
    "TAT": "Y",
    "TAC": "Y",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TAA": "*",
    "TGA": "*",
    "TAG": "*",
}

rev_comp_table = {
    "T": "A",
    "G": "C",
    "C": "G",
    "A": "T",
    "N": "N",
    "Y": "R",
    "R": "Y",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "t": "a",
    "g": "c",
    "c": "g",
    "a": "t",
    "n": "n",
    "y": "r",
    "r": "y",
    "k": "m",
    "m": "k",
    "b": "v",
    "d": "h",
    "h": "d",
    "v": "b",
    "W": "W",
    "S": "S",
    "w": "w",
    "s": "s",
}


def calc_weight(seq):
    seq = seq.replace("*", "")
    weight = aa_weights[seq[0]]
    for char in seq[1:]:
        weight += aa_weights[char] - 18.01
    return weight


def reverse_complement(seq):

    new_seq = ""
    for i in seq[::-1]:
        new_seq += rev_comp_table[i]

    return new_seq


def fq2dict(fastq_input):  # file_ = True):
    fastq_dict = {}
    if fastq_input.endswith((".gz", ".gzip")):
        with gzip.open(fastq_input, "rt") as raw:
            for line in raw:
                data = line.rstrip()
                if data.startswith("@"):
                    header = data[1:].split(" ")
                    seq_name = header[0]
                    fastq_dict[seq_name] = {
                        "sequence": "",
                        "description": " ".join(header[1:]),
                        "score": "",
                    }
                    append_seq = "sequence"
                elif data == "+":
                    append_seq = "score"
                else:
                    fastq_dict[seq_name][append_seq] += data
    else:
        with open(fastq_input, "r") as raw:
            for line in raw:
                data = line.rstrip()
                if data.startswith("@"):
                    header = data[1:].split(" ")
                    seq_name = header[0]
                    fastq_dict[seq_name] = {
                        "sequence": "",
                        "description": " ".join(header[1:]),
                        "score": "",
                    }
                    append_seq = "sequence"
                elif data == "+":
                    append_seq = "score"
                else:
                    fastq_dict[seq_name][append_seq] += data
    return fastq_dict


def dict2fq(fastq_dict, description=True):
    """Convert a fastq dictionary to a fastq string"""
    fastq_string = ""
    if description:
        for seq in fastq_dict:
            fastq_string += ">" + seq.rstrip()
            if "description" in fastq_dict[seq]:
                fastq_string += (" " + fastq_dict[seq]["description"]).rstrip() + "\n"
            else:
                fastq_string += "\n"
            fastq_string += (
                fastq_dict[seq]["sequence"].rstrip()
                + "\n+\n"
                + fastq_dict[seq]["score"].rstrip()
                + "\n"
            )

    else:
        for seq in fastq_dict:
            fastq_string += (
                ">"
                + seq.rstrip()
                + "\n"
                + fastq_dict[seq]["sequence"].rstrip()
                + "\n+\n"
                + fastq_dict[seq]["score"].rstrip()
                + "\n"
            )

    return fastq_string


def xmfa2dict(xmfa_in):
    xmfa_dict = defaultdict(dict)
    count = 0
    with open(xmfa_in, "r") as raw:
        for line in raw:
            data = line.rstrip()
            if data == "=":
                count += 1
            elif data.startswith(">"):
                seq_info = data[2:]
                seq_name, descrip_p = seq_info.split(":")
                coords = " ".join(descrip_p.split()[:2])
                description = " ".join(descrip_p.split()[2:])
                xmfa_dict[seq_name][count] = {
                    "sequence": "",
                    "description": description,
                    "coords": coords,
                }

            elif not data.startswith("#"):
                xmfa_dict[seq_name][count]["sequence"] += data
    if not xmfa_dict[-1]:
        del xmfa_dict[-1]

    return xmfa_dict


def dict2xmfa(xmfa_dict, description=True):
    """Output an xmfa formatted string from xmfa dict"""

    counts = max(list(chain(*[list(x.keys()) for k, x in xmfa_dict.items()])))
    xmfa_string = "#FormatVersion Mauve1\n"
    if description:
        for count in range(counts):
            for gene, genes_dict in xmfa_dict.items():
                if count in genes_dict:
                    gene_dict = genes_dict[count]
                    xmfa_string += f'> {gene.rstrip()}:{gene_dict["coords"]}'
                    if "description" in xmfa_dict[gene]:
                        xmfa_string += f' {gene_dict["description"].rstrip()}\n'
                    else:
                        xmfa_string += "\n"
                    xmfa_string += gene_dict["sequence"].rstrip() + "\n"
            xmfa_string += "=\n"

    return xmfa_string.rstrip()


def fa2dict(fasta_input):  # file_ = True):
    fasta_dict = {}
    with open(fasta_input, "r") as raw:
        for line in raw:
            data = line.rstrip()
            if data.startswith(">"):
                header = data[1:].split(" ")
                seq_name = header[0]
                fasta_dict[seq_name] = {
                    "sequence": "",
                    "description": " ".join(header[1:]),
                }
            elif not data.startswith("#"):
                fasta_dict[seq_name]["sequence"] += data
    return fasta_dict


def fa2dict_str(fasta_input):
    fasta_dict = {}
    for line in fasta_input.split("\n"):
        data = line.rstrip()
        if data.startswith(">"):
            header = data[1:].split(" ")
            seq_name = header[0]
            fasta_dict[seq_name] = {"sequence": "", "description": " ".join(header[1:])}
        elif not data.startswith("#"):
            fasta_dict[seq_name]["sequence"] += data
    return fasta_dict


def fa2dict_accs(fasta_input, accs=set()):  # file_ = True):
    fasta_dict = {}
    with open(fasta_input, "r") as raw:
        for line in raw:
            data = line.rstrip()
            if data.startswith(">"):
                header = data[1:].split(" ")
                seq_name = header[0]
                if seq_name in accs:
                    fasta_dict[seq_name] = {
                        "sequence": "",
                        "description": " ".join(header[1:]),
                    }
                    accs.remove(seq_name)
                elif not accs:
                    break
                else:
                    seq_name = False
            elif seq_name:
                if not data.startswith("#"):
                    fasta_dict[seq_name]["sequence"] += data
    return fasta_dict


# truncates sequences based on inputted lenght
def dnatrunc(fasta_dict, trunc_length):
    for gene in fasta_dict:
        if len(fasta_dict[gene]["sequence"]) >= trunc_length:
            fasta_dict[gene]["sequence"] = fasta_dict[gene]["sequence"][
                : (trunc_length - 1)
            ]
            # I should be able to do this in one line
            fasta_dict[gene]["reverse_complement"] = fasta_dict[gene][
                "reverse_complement"
            ][-1 : (0 - trunc_length) : -1]
            fasta_dict[gene]["reverse_complement"] = fasta_dict[gene][
                "reverse_complement"
            ][::-1]
    return fasta_dict


# takes output dict from fasta2Dict and extracts codons from all 6 reading frames
def dict2codon(fasta_dict):

    fasta_dict[gene]["codons"] = {}

    # outputs a dictionary of reading frames and the possible codons for each
    for index in range(3):
        fasta_dict[gene]["codons"]["reading_frame_" + str(index)] = []
        codon = ""
        for nt in fasta_dict[gene]["sequence"][index:]:
            if len(codon) < 3:
                codon = codon + nt
            if len(codon) == 3:
                fasta_dict[gene]["codons"]["reading_frame_" + str(index)].append(codon)
                codon = ""

    # let's extract codons from all reading frames for the reverse sequences too
    for index in range(3):
        fasta_dict[gene]["codons"]["reverse_reading_frame_" + str(index)] = []
        codon = ""
        for nt in fasta_dict[gene]["reverse_complement"][index:]:
            if len(codon) < 3:
                codon = codon + nt
            if len(codon) == 3:
                fasta_dict[gene]["codons"][
                    "reverse_reading_frame_" + str(index)
                ].append(codon)
                codon = ""

    return codondict


def dict2fa(fasta_dict, description=True):

    fasta_string = ""
    if description:
        for gene in fasta_dict:
            fasta_string += ">" + gene.rstrip()
            if "description" in fasta_dict[gene]:
                fasta_string += (" " + fasta_dict[gene]["description"]).rstrip() + "\n"
            else:
                fasta_string += "\n"
            fasta_string += fasta_dict[gene]["sequence"].rstrip() + "\n"

    else:
        for gene in fasta_dict:
            fasta_string += (
                ">"
                + gene.rstrip()
                + "\n"
                + fasta_dict[gene]["sequence"].rstrip()
                + "\n"
            )

    return fasta_string


# need to turn these into classes and class functions
def calc_gc(gene):

    G = gene["sequence"].count("G")
    C = gene["sequence"].count("C")
    GC = (G + C) / len(gene["sequence"])
    GC_con = "{:.2%}".format(GC)

    return GC_con


class GFFList():
    def __init__():

# GFF strand column -> BioPython-style strand integer ('.'/'?' -> None)
_STRAND = {"+": 1, "-": -1}
# strand integer -> GFF strand column (anything else, e.g. None, serializes to '.')
_STRAND_REVERSE = {1: "+", -1: "-"}
# the placeholder characters GFF3 uses for an undefined numeric/strand column
_UNDEFINED = {".", "?", "", None}

# GFF3 column-9 reserved characters and their percent-encodings; '%' is listed
# first so the escapes introduced below are not themselves re-encoded
_GFF_ATTR_ESCAPES = (
    ("%", "%25"), (";", "%3B"), ("=", "%3D"), ("&", "%26"), (",", "%2C"),
    ("\t", "%09"), ("\n", "%0A"), ("\r", "%0D"),
)

def _gff_escape(value: str) -> str:
    """Percent-encode the GFF3 attribute-reserved characters in ``value``"""
    for char, code in _GFF_ATTR_ESCAPES:
        value = value.replace(char, code)
    return value


def _format_gff_attributes(attributes: dict, field_delimiter: str = ";", value_delimiter: str = "=") -> str:
    """Serialize a ``{key: value}`` attribute dict to a GFF3 column-9 string.

    Reserved characters in keys and values are percent-encoded; an empty dict
    yields ``.`` (the GFF3 "no attributes" placeholder)."""
    if not attributes:
        return "."
    return field_delimiter.join(
        f"{_gff_escape(str(key))}{value_delimiter}{_gff_escape(str(value))}"
        for key, value in attributes.items()
    )


def _as_int(value, name: str):
    """Coerce a GFF numeric column to ``int``, mapping the undefined placeholder
    ('.') to ``None``"""
    if value in _UNDEFINED:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be an integer, got {value!r}")


def _as_strand(value):
    """Coerce a GFF strand column ('+'/'-') or a signed integer to 1/-1, with the
    undefined placeholder ('.'/'?') mapping to ``None``"""
    if value in _UNDEFINED:
        return None
    if value in _STRAND:
        return _STRAND[value]
    if value in (1, -1):
        return value
    raise ValueError(f"strand must be one of '+'/'-'/1/-1, got {value!r}")


class Feature():
    """A single cross-annotation feature
    `start` and `end` are sorted by smallest to largest and are expected to be 0-based half-open upon input"""

    def __init__(self, fid=None, pid=None, seqid=None, source=None, type=None, start=None, end=None, score=None, strand=None, phase=None,
                 attributes: dict = None, descendants: list = None, parent: Feature = None, sequence: str = "", origin_sequence: str = ""
                 ):
        self.id = fid
        # parent ID
        self.pid = pid
        self.seqid = seqid
        self.source = source
        self.type = type
        put_start = _as_int(start, "start")
        put_end = _as_int(end, "end")
        # enforce start < end
        self.start, self.end = sorted((put_start, put_end))
        self.score = _as_int(score, "score")
        self.strand = _as_strand(strand)
        self.phase = _as_int(phase, "phase")
        self.parent = parent

        # sentinel defaults: a shared mutable default would let every Feature
        # alias one dict/list, collapsing the hierarchy built by group_features
        self.attributes = {} if attributes is None else attributes
        self.descendants = [] if descendants is None else descendants

        # derive sequence from the full contiguous sequence provided
        if origin_sequence:
            # 0-BASED, HALF-OPEN
            sequence = origin_sequence[self.start:self.end]

        # demand that the provided sequence abide by the coordinates
        if sequence:
            if len(sequence) < self.end - self.start:
                raise ValueError(f"sequence length deviates from coordinate length")
            else:
                self.sequence = sequence
        else:
            self.sequence = ""


    def _to_gff_line(self) -> str:
        """Serialize this feature (without its descendants) to one GFF3 line.

        Coordinates are converted from the internal 0-based, half-open
        representation back to GFF3's 1-based, both-inclusive columns; an
        undefined ``start``/``end``/``score``/``strand``/``phase`` renders as the
        '.' placeholder."""
        start = "." if self.start is None else str(self.start + 1)
        end = "." if self.end is None else str(self.end)
        score = "." if self.score is None else str(self.score)
        strand = _STRAND_REVERSE.get(self.strand, ".")
        phase = "." if self.phase is None else str(self.phase)
        id_dict = {"ID": self.fid}
        if self.pid:
            id_dict["Parent"] = self.pid
        return "\t".join((
            self.seqid, self.source, self.type, start, end,
            score, strand, phase, _format_gff_attributes({**id_dict, **self.attributes)),
        ))
    
    def to_gff(self) -> str:
        """Serialize this feature and all of its descendants to a GFF3 string.

        Lines are emitted depth-first, each parent before its children, joined
        by newlines (no trailing newline)."""
        lines = [self._to_gff_line()]
        for descendant in self.descendants:
            lines.append(descendant.to_gff())
        return "\n".join(lines)


def _attr_get(attributes: dict, keys):
    """Return the value of the first present attribute key (accepts case
    variants, e.g. ``ID``/``id`` or ``Parent``/``parent``), or ``None``"""
    for key in keys:
        if key in attributes:
            return attributes[key]
    return None


def group_features(features: list, parent_ids: list = ["Parent", "parent"], ids: list = ["ID", "Id", "id"]):
    """Hierarchically group features based on Parent <-> ID relationships"""
    id2feature = {}
    for feature in features:
        fid = _attr_get(feature.attributes, ids)
        if fid in id2feature:
            raise KeyError(f"{fid} is depicted in multiple Features")
        id2feature[fid] = feature

    feature_dict = defaultdict(list)
    for fid, feature in id2feature.items():
        par_id = _attr_get(feature.attributes, parent_ids)
        if par_id:
            # link features together
            id2feature[par_id].descendants.append(feature)
            feature.parent = id2feature[par_id]
        feature_dict[seqid].append(feature)

    return feature_dict






# need to change into a class
def gff2list(gff_info, path=True, error=True):

    gff_list_dict = []
    data = []
    if path:
        with open(gff_info, "r") as raw_gff:
            for line in raw_gff:
                d = line.rstrip()
                if d and not d.startswith("#"):
                    data.append(d.split("\t"))
                elif d.startswith("##FASTA"):
                    break
    else:
        for line in gff_info.split("\n"):
            if not line.startswith("#") and line:
                d = line.split("\t")
                data.append(d)
            elif line.startswith("##FASTA"):
                break
    try:
        for col_list in data:
            gff_list_dict.append(
                {
                    "seqid": col_list[0],
                    "source": col_list[1],
                    "type": col_list[2],
                    "start": int(col_list[3]),
                    "end": int(col_list[4]),
                    "score": col_list[5],
                    "strand": col_list[6],
                    "phase": col_list[7],
                    "attributes": col_list[8],
                }
            )
    except IndexError:
        raise IndexError(
            str(len(col_list))
            + "/9 expected tab-"
            + "delimitted fields: "
            + str(col_list)
        )
    except ValueError:
        raise ValueError(
            str(col_list[4:6]) + " invalid integer " + "conversion: " + str(col_list)
        )

    return gff_list_dict


def list2gff(gff_list, ver=3):

    if ver:
        gff_str = "##gff-version " + str(ver) + "\n"
    else:
        gff_str = ""
    for line in gff_list:
        add_str = "\t".join(str(val) for val in list(line.values()))
        gff_str += add_str + "\n"

    return gff_str.rstrip()


def gff3_comps(source=None):

    comps = {}
    comps["par"] = "(?:^|(?<=;))" + r'Parent=["\']?([^;\'"]+)'
    comps["id"] = "(?:^|(?<=;))" + r'ID=["\']?([^;\'"]+)'
    comps["Alias"] = "(?:^|(?<=;))" + r'Alias=["\']?([^;\'"]+)'
    comps["product"] = "(?:^|(?<=;))" + r"""product=["\']?([^;"']+)"""
    comps["OG"] = "(?:^|(?<=;))" + r'OG=["\']?([\w+:\d+\|]+)'
    comps["ver"] = "gff3"
    comps["prot"] = (
        "(?:^|(?<=;))"
        + r'protein_id=["\']?([^;\'"]+)|'
        + "(?:^|(?<=;))"
        + r'proteinId=["\']?([^"\';]+)'
    )
    comps["transcript"] = (
        "(?:^|(?<=;))"
        + r'transcriptId=["\']?([^;"\'])|'
        + "(?:^|(?<=;))"
        + r'transcript_id=["\']?([^;\'"]+)'
    )

    return comps


def gff2_comps():

    comps = {}
    comps["id"] = r'name "([^"]+)"'
    # this unfortunately includes "gene_name" from gtf naming
    comps["prot"] = r"proteinId ([^;]+)"
    comps["transcript"] = r"transcriptId ([^;]+)"
    comps["alias"] = r'alias "([^"]+)"'
    comps["product"] = r'product_name "([^"]+)'
    comps["ver"] = "gff2"

    return comps


def gtf_comps():

    comps = {}
    comps["id"] = r'gene_id "?([^"]+)"?'
    comps["transcript"] = r'transcript_id "([^"]+)"'
    comps["alias"] = r'alias "([^"]+)"'
    comps["ver"] = "gtf"

    return comps


def compile_exon(gff):

    exon_dict = {}

    for index in range(len(gff)):
        if gff[index]["type"].lower() == "exon":
            break

    protComp = re.compile(r";Parent\=([^;]*)")
    if not protComp.search(gff[index]["attributes"]):
        protComp = re.compile(r'gene_id "(.*?)"')
        if not protComp.search(gff[index]["attributes"]):
            protComp = re.compile(r'name "(.*?)"\;')
            if not protComp.search(gff[index]["attributes"]):
                protComp = re.compile(r"ID=(.*?);")

    for line in gff:
        if line["type"].lower() == "exon":
            prot = protComp.search(line["attributes"])[1]
            if prot not in exon_dict:
                exon_dict[prot] = []
            if line["strand"] == "+":
                exon_dict[prot].append([int(line["start"]) - 1, int(line["end"])])
            else:
                if int(line["start"]) > int(line["end"]):
                    exon_dict[prot].append([int(line["end"]) - 1, int(line["start"])])
                else:
                    exon_dict[prot].append([int(line["start"]) - 1, int(line["end"])])

    for prot in exon_dict:
        exon_dict[prot] = sorted(exon_dict[prot], key=lambda i: i[0])

    return exon_dict

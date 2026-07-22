#! /usr/bin/env python3
"""Offline tests for JGI Data Portal file selection - `download jgi`.

Focused on keeping mitochondrial assemblies out of the organismal genome slot.
MycoCosm files its mitochondrial assemblies under a *nuclear* label: they carry
``jat_label`` ``assembly_unmasked`` and sit in "Genome Assembly (unmasked)",
exactly like the real genome (a survey of 250 portals found 93 such files, and
Lst7536_1 lists one file under both the nuclear and the ``assembly_mitochondrial``
label). So the label cannot discriminate and `select_file` screens the filename.

The substitution is silent and destructive - a mitochondrion is ~1/1000th the
size of the genome it would replace, so it yields a nonsense MTDB entry rather
than an obvious failure. Every record below is copied from a live JGI listing.
"""
import pytest

from mycotools.download.jgi import _is_mito, select_file


def _f(name, label, fmt="fasta", status="RESTORED"):
    """Build the file record shape `select_file` consumes."""
    return {
        "file_name": name,
        "file_status": status,
        "metadata": {"jat_label": label, "file_format": fmt},
    }


# --------------------------------------------------------------------------- #
# _is_mito - filename screen
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "name",
    [
        # mislabelled as nuclear (assembly_unmasked) - the actual bug
        "Suilu4_MitoAssemblyScaffolds.fasta.gz",
        "Amnli1_MitoScaffolds.fasta.gz",
        "Lst7536_1_MitoAssemblyScaffolds_2023-09-19.fasta.gz",
        "Umbsp_AD052_1_MitoAssemblyScaffolds.fasta.gz",
        # correctly labelled assembly_mitochondrial - screened anyway
        "Agrped1_MitochondrialScaffolds.fasta.gz",
        "Copmic2_MitochondrionScaffolds.fasta.gz",
        "Xxx1_MitochondriumScaffolds.fasta.gz",
        "Xxx1_MitochondrianScaffold.fasta",
        "Xxx1_mitochondrial_scaffold.fasta.gz",
        # mitochondrial annotations, labelled genes_filtered
        "Lst7536_1_MitoGenes_2023-09-19.gff3.gz",
        # separator and position variants
        "Mito_assembly.fasta.gz",
        "Xxx1-MitoContigs.fa",
    ],
)
def test_mitochondrial_files_are_recognized(name):
    assert _is_mito({"file_name": name})


@pytest.mark.parametrize(
    "name",
    [
        # real nuclear assemblies, across MycoCosm's naming eras
        "Neucr2_AssemblyScaffolds_Repeatmasked.fasta.gz",
        "Suilu4_AssemblyScaffolds.fasta.gz",
        "Cenge1058_1_AssemblyScaffolds_Repeatmasked_2023-10-07.fasta.gz",
        "Heterobasidion_annosum.AssembledScaffolds.fasta.gz",
        "Abisporus_varbisporusH97.v2_AssembledScaffolds.fasta.gz",
        "PleosPC9_1_assembly_scaffolds_repeatmasked.fasta.gz",
        "Neucr2_GeneCatalog_genes_20130412.gff.gz",
        # organism names that merely *contain* "mito" - the false-positive trap.
        # Fomitopsis carries it mid-token; Mitosporidium is a real fungal genus
        # that carries it token-initially, and full organism names do reach
        # filenames (see the Heterobasidion entry above).
        "Fomitopsis_schrenkii_FP58527_AssemblyScaffolds.fasta.gz",
        "Fomitiporia_mediterranea_AssemblyScaffolds.fasta.gz",
        "Mitosporidium_daphniae_AssemblyScaffolds.fasta.gz",
        "Fompi3_AssemblyScaffolds.fasta.gz",
    ],
)
def test_organismal_files_are_not_mistaken_for_mitochondrial(name):
    assert not _is_mito({"file_name": name})


# --------------------------------------------------------------------------- #
# select_file - the mitochondrion must never fill the genome slot
# --------------------------------------------------------------------------- #
# Suilu4 as JGI lists it: the mitochondrion is RESTORED and shares the
# unmasked label with the real genome, so `download jgi -a -n` picked the 33 KB
# mitochondrion over the 13.7 MB genome.
SUILU4 = [
    _f("Suilu4_MitoAssemblyScaffolds.fasta.gz", "assembly_unmasked"),
    _f("Suilu4_AssemblyScaffolds_Repeatmasked.fasta.gz", "assembly_masked"),
    _f("Suilu4_AssemblyScaffolds.fasta.gz", "assembly_unmasked"),
]


@pytest.mark.parametrize("masked", [True, False], ids=["masked", "nonmasked"])
def test_genome_is_chosen_over_mitochondrion(masked):
    chosen = select_file(SUILU4, "fna", masked=masked)
    expect = (
        "Suilu4_AssemblyScaffolds_Repeatmasked.fasta.gz"
        if masked
        else "Suilu4_AssemblyScaffolds.fasta.gz"
    )
    assert chosen["file_name"] == expect


@pytest.mark.parametrize("masked", [True, False], ids=["masked", "nonmasked"])
def test_choice_does_not_depend_on_jgi_listing_order(masked):
    """The mitochondrion and the unmasked genome share a label, so before the
    filename screen the winner was decided by a stable-sort tie - i.e. by the
    order JGI happened to return. `mtdb update` (masked=True) was one
    reordering away from storing mitochondria."""
    for order in (SUILU4, list(reversed(SUILU4)), SUILU4[1:] + SUILU4[:1]):
        assert not _is_mito(select_file(order, "fna", masked=masked))


def test_tape_status_cannot_promote_a_mitochondrion():
    """RESTORED-over-PURGED outranks label preference, and mitochondria are
    small enough to stay on disk while the genome is archived. That must not
    be enough to select one."""
    files = [
        _f("Boled5_MitoAssemblyScaffolds.fasta.gz", "assembly_unmasked"),
        _f("Boled5_AssemblyScaffolds.fasta.gz", "assembly_unmasked", status="PURGED"),
        _f(
            "Boled5_AssemblyScaffolds_Repeatmasked.fasta.gz",
            "assembly_masked",
            status="PURGED",
        ),
    ]
    chosen = select_file(files, "fna", masked=True)
    assert chosen["file_name"] == "Boled5_AssemblyScaffolds_Repeatmasked.fasta.gz"


def test_mitochondrial_annotations_are_excluded_from_gff3():
    """`*_MitoGenes_*.gff3.gz` carries genes_filtered - the same mislabelling,
    and a mitochondrial gff3 against a nuclear assembly is a corrupt entry."""
    files = [
        _f("Lst7536_1_MitoGenes_2023-09-19.gff3.gz", "genes_filtered", fmt="gff3"),
        _f(
            "Lst7536_1_FilteredModels1_2023-09-10.gff3.gz", "genes_filtered", fmt="gff3"
        ),
    ]
    for order in (files, list(reversed(files))):
        chosen = select_file(order, "gff3")
        assert chosen["file_name"] == "Lst7536_1_FilteredModels1_2023-09-10.gff3.gz"


def test_mitochondrion_only_organism_yields_no_assembly():
    """Returning None marks the organism failed; it must not enter MTDB with a
    mitochondrion standing in for its genome."""
    files = [_f("Xxx1_MitoAssemblyScaffolds.fasta.gz", "assembly_unmasked")]
    assert select_file(files, "fna", masked=True) is None
    assert select_file(files, "fna", masked=False) is None


# --------------------------------------------------------------------------- #
# regression guard - ordinary organisms are untouched
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "masked,expect",
    [
        (True, "Neucr2_AssemblyScaffolds_Repeatmasked.fasta.gz"),
        (False, "Neucr2_AssemblyScaffolds.fasta.gz"),
    ],
)
def test_masked_preference_is_preserved(masked, expect):
    files = [
        _f("Neucr2_AssemblyScaffolds_Repeatmasked.fasta.gz", "assembly_masked"),
        _f("Neucr2_AssemblyScaffolds.fasta.gz", "assembly_unmasked"),
    ]
    assert select_file(files, "fna", masked=masked)["file_name"] == expect


def test_restored_still_outranks_purged_among_genomes():
    """The tape-avoidance preference survives the new screen."""
    files = [
        _f(
            "Neucr2_AssemblyScaffolds_Repeatmasked.fasta.gz",
            "assembly_masked",
            status="PURGED",
        ),
        _f("Neucr2_AssemblyScaffolds.fasta.gz", "assembly_unmasked"),
    ]
    chosen = select_file(files, "fna", masked=True)
    assert chosen["file_name"] == "Neucr2_AssemblyScaffolds.fasta.gz"

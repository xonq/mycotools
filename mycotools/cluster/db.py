#! /usr/bin/env python3

import sys
import shutil
import logging
import argparse
import subprocess
import multiprocessing as mp
from statistics import stdev, StatisticsError
from collections import defaultdict, Counter
from mycotools.mtdb.files import soft_main as symlink_files
from mycotools.lib.dbtools import mtdb, primary_db
from mycotools.lib.biotools import dict2fa, fa2dict_accs
from mycotools.lib.kontools import format_path, mk_output, find_execs, setup_logging
from pathlib import Path


logger = logging.getLogger(__name__)


def mk_db2hg_output(out_dir, nscg=False):
    """Prepare the output directories"""
    wrk_dir = out_dir + "working/"
    scg_dir = out_dir + "single_copy_genes/"
    nscg_dir = out_dir + "near_single_copy_genes/"
    hg_seq_dir = out_dir + "hgs/"

    if not Path(wrk_dir).is_dir():
        Path(wrk_dir).mkdir()
    if not Path(scg_dir).is_dir():
        Path(scg_dir).mkdir()
    if not Path(nscg_dir).is_dir() and nscg:
        Path(nscg_dir).mkdir()
    if not Path(hg_seq_dir).is_dir():
        Path(hg_seq_dir).mkdir()

    return wrk_dir, scg_dir, nscg_dir, hg_seq_dir


def run_mmseqs(
    db,
    wrk_dir,
    algorithm="mmseqs easy-cluster",
    min_id=0.3,
    min_cov=0.3,
    sensitivity=7.5,
    cpus=1,
):
    """Run MMseqs clustering by sym linking MTDB proteomes"""
    symlink_files(["faa"], db, wrk_dir, verbose=False)  # symlink proteomes
    cluster_res_file = wrk_dir + "raw_hgs.tsv"
    if not Path(cluster_res_file).is_file():  # NEED to add to log removal
        # be cautious about shell injection because we need to glob
        int(cpus)
        float(min_id)
        float(min_cov)
        if not Path(wrk_dir).is_dir():
            raise OSError("invalid working directory")
        elif not algorithm in {"mmseqs easy-linclust", "mmseqs easy-cluster"}:
            raise OSError("invalid mmseqs binary")
        mmseqs_cmd = subprocess.call(
            algorithm
            + " "
            + " ".join(
                [
                    wrk_dir + "faa/*faa",
                    wrk_dir + "cluster",
                    wrk_dir + "tmp/",
                    "--min-seq-id",
                    str(min_id),
                    "--threads",
                    str(cpus),
                    "--compressed",
                    "1",
                    "-s",
                    str(sensitivity),
                    "--cov-mode",
                    "0",
                    "-c",
                    str(min_cov),
                ]
            ),
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        #                                      stderr = subprocess.DEVNULL)
        shutil.move(wrk_dir + "cluster_cluster.tsv", cluster_res_file)
    elif Path(cluster_res_file).stat().st_size:
        mmseqs_cmd = 0
    else:
        mmseqs_cmd = 1
    if mmseqs_cmd:
        logger.error("cluster failed")
        sys.exit(1)
    if Path(wrk_dir + "cluster_all_seqs.fasta").is_file():
        Path(wrk_dir + "cluster_all_seqs.fasta").unlink()
    if Path(wrk_dir + "cluster_rep_seq.fasta").is_file():
        Path(wrk_dir + "cluster_rep_seq.fasta").unlink()
    if Path(wrk_dir + "tmp/").is_dir():
        shutil.rmtree(wrk_dir + "tmp/")
    return cluster_res_file


def parse_1to1(hg_file, useableOmes=set()):
    """Parse MMseqs style homology group output"""
    derivations = defaultdict(list)
    with open(hg_file, "r") as raw:
        for line in raw:
            k, v = line.rstrip().split()  # all white space
            derivations[k].append(v)

    hgs_list = []
    for k, v in derivations.items():
        v.append(k)
        hgs_list.append(sorted(set(v)))

    hg2gene = {
        i: v for i, v in enumerate(sorted(hgs_list, key=lambda x: len(x), reverse=True))
    }

    for hg, genes in hg2gene.items():
        genes = [x for x in genes if x[: x.find("_")] in useableOmes]
    hg2gene = {
        k: v for k, v in sorted(hg2gene.items(), key=lambda x: len(x[1]), reverse=True)
    }

    gene2hg = {}
    for i, hg_list in hg2gene.items():
        for gene in hg_list:
            gene2hg[gene] = i

    i2ome = sorted(set([x[: x.find("_")] for x in list(gene2hg.keys())]))
    ome_num = {v: i for i, v in enumerate(i2ome)}

    return ome_num, gene2hg, i2ome, hg2gene


def compile_homolog_groups(raw_hg_file, out_hg_file, useableOmes=set()):
    """Compile OrthoFinder-style homology group tabular output"""
    hg_info = parse_1to1(raw_hg_file)
    hg2genes = hg_info[-1]
    with open(out_hg_file, "w") as out:
        for hg, genes in hg2genes.items():
            out.write(str(hg) + "\t" + " ".join(genes) + "\n")

    return hg_info


def id_near_schgs(
    hg2gene, omes, max_hgs=10000, max_median=100, max_mean=2, max_stdev=1, min_genomes=1
):
    """Identify near single copy homology groups"""
    schgs, near_schgs = [], []
    hg2stats, all_omes_hgs, hg2d_omes = {}, set(), {}
    min_hg2gene = {
        k: v
        for k, v in sorted(hg2gene.items(), key=lambda x: len(x[1]))
        if len(v) >= len(omes)
    }

    if max_mean:
        max_len = max_mean * len(omes)  # dont need to compute means iteratively
    else:
        max_len = None

    for hg, genes in min_hg2gene.items():
        hg_omes = [x[: x.find("_")] for x in genes]
        count_omes = list(Counter(hg_omes).values())  # prepare for median
        try:
            sd = stdev(count_omes)
        except StatisticsError:
            sd = None
        # calculation
        count_omes.sort()
        median_i = round((len(count_omes) / 2) - 0.5)
        median = count_omes[median_i]
        set_hg_omes = set(hg_omes)
        omes_perc = len(set_hg_omes) / len(omes)
        hg2stats[hg] = (
            len(genes) / len(omes),
            median,
            sd,
            omes_perc,
        )
        diff_omes = omes.difference(set_hg_omes)
        hg2d_omes[hg] = diff_omes

        if omes_perc >= min_genomes:  # all omes are present
            if omes_perc == 1:
                all_omes_hgs.add(hg)
            #            print(median, max_median, len(genes), max_len)
            if max(count_omes) == 1:
                schgs.append(hg)
                near_schgs.append(hg)
            else:
                if max_median:
                    if median > max_median:
                        continue
                if max_len:
                    if len(genes) > max_len:
                        break  # it is sorted so nothing will be lower
                if max_stdev:
                    if sd is None:
                        continue
                    elif sd > max_stdev:
                        continue
                near_schgs.append(hg)
    #            if len(near_schgs) == max_hgs:
    #               break
    return schgs, near_schgs, hg2stats, all_omes_hgs, hg2d_omes


def write_hg_stats(hg2stats, out_file, sort=False):
    """Write the mean, median, and stdev for the stats"""
    if sort:
        # sort by median, then standard dev
        hg2stats = {
            k: v
            for k, v in sorted(
                hg2stats.items(), key=lambda x: (x[1][1], (x[1][2] is None, x[1][2]))
            )
        }

    with open(out_file, "w") as out:
        out.write("#hg\tmean\tmedian\tstdev\t%_genomes\n")
        for hg, stats in hg2stats.items():
            out.write(f"{hg}\t{stats[0]}\t{stats[1]}\t{stats[2]}" + f"\t{stats[3]}\n")


def write_hgs(hg, genes, wrk_dir, write_dir):
    """Write homology group fastas"""
    #    fa_dict = acc2fa(db, genes)
    ome2gene = defaultdict(list)
    for gene in genes:
        ome = gene[: gene.find("_")]
        ome2gene[ome].append(gene)

    fa_dict = {}
    for ome, genes in ome2gene.items():
        fa_dict = {**fa_dict, **fa2dict_accs(f"{wrk_dir}faa/{ome}.faa", set(genes))}

    with open(f"{write_dir}{hg}.faa.tmp", "w") as out:
        out.write(dict2fa(fa_dict))
    Path(f"{write_dir}{hg}.faa.tmp").rename(f"{write_dir}{hg}.faa")


def align_hg(hg_fa, out_fa, cpus=1):
    """Run homology group alignment"""
    with open(out_fa + ".tmp", "w") as out:
        cmd = subprocess.call(
            ["mafft", "--auto", "--thread", f"-{cpus}", hg_fa], stdout=out
        )  # , stderr = subprocess.PIPE)
    Path(out_fa + ".tmp").rename(out_fa)
    return cmd


def hmmbuild_hg(msa_fa, out_hmm, cpus=1):
    """Build HMMs of homology groups"""
    cmd = subprocess.call(
        ["hmmbuild", "--amino", out_hmm + ".tmp", msa_fa],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    Path(out_hmm + ".tmp").rename(out_hmm)
    return cmd


def pangenome_output(pan_file, aln_file, hg2gene, hg2d_omes, max_mis_ome=0):
    """Circumscribe and output the pangenome based on a maximum missing genome
    to determine the accessory and core HGs"""
    ome2pan = defaultdict(lambda: [[], []])
    core_hgs = set(k for k, v in hg2d_omes.items() if len(v) <= max_mis_ome)
    #    acc_hgs = sorted(set(hg2d_omes.keys()).difference(core_hgs))
    acc_hgs = set()
    ome2acc_aln = defaultdict(set)
    for hg, genes in hg2gene.items():
        if hg in core_hgs:
            for gene in genes:
                ome = gene[: gene.find("_")]
                ome2pan[ome][0].append(gene)
        else:
            acc_hgs.add(hg)
            for gene in genes:
                ome = gene[: gene.find("_")]
                ome2pan[ome][1].append(gene)
                ome2acc_aln[ome].add(hg)

    acc_hgs = sorted(acc_hgs)
    ome2aln = {}
    for ome, hg_set in ome2acc_aln.items():
        ome2aln[ome] = {k: 0 for k in acc_hgs}
        for hg in list(hg_set):
            ome2aln[ome][hg] = 1

    with open(pan_file, "w") as out:
        out.write("#ome\tcore_%\tacc_%\n")
        for ome, pan in ome2pan.items():
            tot = len(pan[0]) + len(pan[1])
            out.write(f"{ome}\t{len(pan[0])/tot}\t{len(pan[1])/tot}\n")

    with open(aln_file, "w") as out:
        #        out.write('#ome ' + ' '.join([str(x) for x in acc_hgs]) + '\n')
        out.write(f"{len(ome2aln)} {len(acc_hgs)}\n")
        for ome, aln in {
            k: v for k, v in sorted(ome2aln.items(), key=lambda x: x[0])
        }.items():
            out.write(f"{ome} " + "".join([str(x) for x in aln.values()]) + "\n")

    return ome2pan


def main(
    db,
    out_dir,
    min_id=0.3,
    min_cov=0.3,
    sensitivity=7.5,
    max_mean=1,
    max_stdev=1,
    all_hgs=False,
    min_genomes=1,
    max_median=1,
    algorithm="mmseqs easy-cluster",
    hmm=False,
    cpus=1,
):
    """Acquire the proteome files for an inputted MycotoolsDB, use MMseqs to
    cluster protein sequences based on sequence similarity, and output near
    single-copy gene groups"""

    db = db.set_index()

    nscg = False
    if max_mean:
        if max_mean > 1:
            nscg = True
    if max_median and not nscg:
        if max_median > 1:
            nscg = True
    if max_stdev and not nscg:
        if max_stdev > 0:
            nscg = True
    if min_genomes and not nscg:
        if min_genomes < 1:
            nscg = True

    hg_output = out_dir + "homology_groups.tsv"
    hg_stats_file = out_dir + "hg_stats.tsv"
    full_hg_stats_file = out_dir + "core_hg_stats.tsv"
    hg2missing_genome_file = out_dir + "hg2missing.tsv"
    pan_file = out_dir + "pangenome.tsv"
    aln_file = out_dir + "accessory_alignment.phy"
    wrk_dir, scg_dir, nscg_dir, hg_seq_dir = mk_db2hg_output(out_dir, nscg)

    logger.info("Clustering protein sequences")
    raw_hg_file = run_mmseqs(db, wrk_dir, algorithm, min_id, min_cov, sensitivity, cpus)

    useable_omes = set(db.keys())
    logger.info("Compiling homology groups (HGs)")
    ome_num, gene2hg, i2ome, hg2gene = compile_homolog_groups(
        raw_hg_file, hg_output, useable_omes
    )

    logger.info("Identifying single-copy HGs")
    schgs, nschgs, hg2stats, full_hgs, hg2d_omes = id_near_schgs(
        hg2gene,
        set(i2ome),
        max_hgs=100,
        max_median=max_median,
        max_mean=max_mean,
        max_stdev=max_stdev,
        min_genomes=min_genomes,
    )

    logger.info("Writing output")
    pangenome_output(pan_file, aln_file, hg2gene, hg2d_omes, max_mis_ome=0)

    with open(hg2missing_genome_file, "w") as out:
        out.write("#hg\tmissing\n")
        hg2d_omes = {
            k: v for k, v in sorted(hg2d_omes.items(), key=lambda x: len(x[1]))
        }
        for k, v in hg2d_omes.items():
            out.write(f'{k}\t{",".join(sorted(v))}\n')
    write_hg_stats(hg2stats, hg_stats_file, sort=True)
    write_hg_stats(
        {k: hg2stats[k] for k in list(full_hgs)}, full_hg_stats_file, sort=True
    )
    if nscg:
        logger.info(f"{len(nschgs)} near single-copy HGs")
        with mp.Pool(processes=cpus) as pool:
            pool.starmap(
                write_hgs,
                (
                    (hg, hg2gene[hg], wrk_dir, nscg_dir)
                    for hg in nschgs
                    if not Path(f"{nscg_dir}{hg}.faa").is_file()
                ),
            )
    #        for hg in nschgs:
    #           if not os.path.isfile(f'{nscg_dir}{hg}.faa'):
    #              write_hgs(db, hg, hg2gene[hg], wrk_dir, nscg_dir)
    if schgs:
        logger.info(f"{len(schgs)} single-copy HGs")
        with mp.Pool(processes=cpus) as pool:
            pool.starmap(
                write_hgs,
                (
                    (hg, hg2gene[hg], wrk_dir, scg_dir)
                    for hg in schgs
                    if not Path(f"{scg_dir}{hg}.faa").is_file()
                ),
            )
    #        for hg in schgs:
    #           if not os.path.isfile(f'{scg_dir}{hg}.faa'):
    #              write_hgs(db, hg, hg2gene[hg], wrk_dir, scg_dir)
    if all_hgs:
        with mp.Pool(processes=cpus) as pool:
            pool.starmap(
                write_hgs,
                (
                    (hg, genes, wrk_dir, hg_seq_dir)
                    for hg, genes in hg2genes.items()
                    if not Path(f"{hg_seq_dir}{hg}.faa").is_file()
                ),
            )

    #        for hg, genes in hg2gene.items():
    #           if not os.path.isfile(f'{hg_seq_dir}{hg}.faa'):
    #              write_hgs(db, hg, genes, wrk_dir, hg_seq_dir)

    if hmm:
        logger.info("Aligning and building HMMs")
        msa_dir = out_dir + "msa/"
        hmm_dir = out_dir + "hmm/"
        for d in [msa_dir, hmm_dir]:
            if not Path(d).is_dir():
                Path(d).mkdir()
        if nscg:
            srch_hgs = nschgs
        else:
            srch_hgs = schgs
        for hg in srch_hgs:
            if not Path(f"{msa_dir}{hg}.mafft.faa").is_file():
                align_hg(
                    f"{nscg_dir}{hg}.faa", f"{msa_dir}{hg}.mafft.faa", cpus=cpus
                )
            if not Path(f"{hmm_dir}{hg}.hmm").is_file():
                hmmbuild_hg(
                    f"{msa_dir}{hg}.mafft.faa", f"{hmm_dir}{hg}.hmm", cpus=cpus
                )


def cli():
    parser = argparse.ArgumentParser(
        description="Circumscribe protein sequences into homology groups"
    )
    parser.add_argument("-d", "--mtdb", default=primary_db())
    parser.add_argument(
        "-m",
        "--mean",
        type=float,
        help="Maximum mean genes/HG for near single-copy HGs",
    )
    parser.add_argument(
        "-e",
        "--median",
        type=float,
        help="Maximum median genes/HG for near single-copy HGs",
    )
    parser.add_argument(
        "-s",
        "--stdev",
        type=float,
        help="Maximum stdev genes/HG for near single-copy HGs",
    )
    parser.add_argument(
        "-g",
        "--min_genomes",
        type=float,
        help="Minimum percent genomes for near single-copy HGs " + "0 < -g < 1",
        default=1,
    )
    parser.add_argument(
        "-a", "--all", action="store_true", help="Write all HG sequences"
    )
    parser.add_argument(
        "--hmm",
        action="store_true",
        help="Align & build HMMER models of near single-copy HGs",
    )
    parser.add_argument("-o", "--out_dir", help="WARNING: will not overwrite")
    parser.add_argument("-c", "--cpus", type=int, default=mp.cpu_count())
    args = parser.parse_args()
    setup_logging(verbose=getattr(args, "verbose", False))

    if args.min_genomes > 1 or args.min_genomes < 0:
        logger.error("--min_genomes must be between 0 and 1")
        sys.exit(143)

    db = mtdb(format_path(args.mtdb))

    if args.out_dir:
        out_dir = format_path(args.out_dir)
        if not Path(out_dir).is_dir():
            Path(out_dir).mkdir()
            out_dir += "/"
    else:
        out_dir = mk_output(args.out_dir, "cluster_db")

    execs = ["mmseqs"]
    if args.hmm:
        execs.extend(["hmmbuild", "mafft"])
    find_execs(execs, set(execs))

    main(
        db,
        out_dir,
        min_id=0.3,
        min_cov=0.3,
        sensitivity=7.5,
        max_mean=args.mean,
        max_stdev=args.stdev,
        all_hgs=args.all,
        algorithm="mmseqs easy-cluster",
        min_genomes=args.min_genomes,
        hmm=args.hmm,
        max_median=args.median,
        cpus=args.cpus,
    )


if __name__ == "__main__":
    cli()
    sys.exit(0)

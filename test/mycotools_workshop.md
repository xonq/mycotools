## Prerequisites

- A UNIX/Linux command-line environment. Mycotools is deveoped for Linux,
  though it should be compatible with UNIX-based BASH environments, Linux
  virtual machines, and potentially Windows Subsystem for Linux. UNFORTUNATELY, MACS
  THAT USE THE NEW M* SERIES CPUS DO NOT HAVE NATIVE PACKAGE SUPPORT FOR MANY OF THE
  DEPENDENCIES. If you are using a Mac that has an Intel CPU, you will need to 
  switch to BASH as your default shell by opening a terminal and running: `chsh -s /bin/bash`
- [Install Mycotools](https://github.com/xonq/mycotools/blob/master/README.md#install) using miniconda3 as an environment manager
- [Install FigTree](https://github.com/rambaut/figtree/releases) for
  phylogenetic tree viewing
- Basic bash shell navigation knowledge
- [Setup an NCBI account](https://account.ncbi.nlm.nih.gov/)
- [Setup a JGI account](https://signon.jgi.doe.gov/signon) and review the
  use-restricted data [terms and conditions](https://jgi.doe.gov/user-programs/pmo-overview/policies/legacy-data-policies/)

<br />

### Optional

- Install Clinker in your activated Mycotools environment (`conda activate
  mycotools`) via `conda install clinker-py -c conda-forge -c bioconda`
- Install OrthoFinder in your activated Mycotools environment via `conda
  install orthofinder -c bioconda`

<br /><br />

## Make a directory for our analyses
One of the best habits while programming is to make clean,
organized directory hierarchies for projects. I personally recommend starting
an overarching `projects` directory, and any subdirectories titled
`<PROJECT_NAME>_YYYYmm` where `YYYYmm` is the year and month. Navigate
somewhere you want to start your directories (Desktop or Home) and initialize
the hierarchy:

```bash
mkdir <PROJECT_DIR>/mycotools_202405
cd <PROJECT_DIR>/mycotools_202405
```

<br />

## Download the reference database
```bash
curl -o reference.mtdb https://raw.githubusercontent.com/xonq/mycotools/master/test/reference.mtdb
```

<br />

## Initializing a MycotoolsDB (MTDB) from a reference dataset (.mtdb file)
There are multiple ways to initialize a Mycotools database (MTDB): from a de
novo assimilation of publicly available data (`mtdb update -i <INIT_DIR>`);
from a spreadsheet of metadata referencing local genomes (`mtdb update -i
<INIT_DIR> --predb <PREDB.tsv>`); or from a reference MTDB file that contains
publicly available assembly accessions. We will implement a reference.

```bash
# access mycotools by making sure your environment is activated
conda activate mycotools

# initialize the mycotoolsdb according to a reference file
mtdb update -i <PROJECT_DIR>/mycotools_202405/ -r reference.mtdb
```

<br />

## Add a local genome to this database
If you have your own genomes, we can add those to the database by filling out a
spreadsheet of metadata that references those genomes. In this example, we will
add a local genome that we've acquired from JGI - note that Mycotools can
automatically assimilate JGI genomes by initializing the database with `mtdb
update -i <INIT_DIR>`, but that will take quite some time. So we will take one genome for this
example. We will need both an assembly and a `gff3` annotation file to add a
genome to the database.

```bash
# make a directory to run these commands
mkdir jgi_dwnld_202405
cd jgi_dwnld_202405

# download the assembly and gff of a MycoCosm accession
jgiDwnld -i Ustbr1 -a -g
```

This script will output a `predb` file that is ready for assimilating into the
database. If you wanted to add your own genomes, you would fill out one of
these files manually by generating a blank copy via `mtdb predb2mtdb > predb.tsv`,
then running the following commands as we will here:

```bash
# curate the data via predb2mtdb
mtdb p Ustbr1.predb.tsv

# add the curated data to the primary MTDB
mtdb u -a predb2mtdb_<YYYYmmdd>/predb2mtdb.mtdb
```

Now we can check if the file was added by querying the genome code from the
database:

```bash
mtdb ustbro1
```

<br />

## A note on the `mtdb` database command
`mtdb` is a command that controls your interface to generated MTDBs. You can have
different *primary* MTDBs for different projects by initializing different
databases, and/or you can have a large primary database and interact with subsets
of it/

To print the path of the primary database you are linked to, simply run:

```bash
mtdb
```

That's it! Now, we can interface with that primary data using some basic shell
wizardy, e.g. to preview the contents of the primary MTDB:

```bash
cat $(mtdb)
```

Now, if you want to interface with a subset of your primary MTDB we can
*extract* the particular genomes of interest:

```bash
# extract a database of Ustilago genomes
mtdb extract -l Ustilago > ust.mtdb
```

<br />

## Obtain basic genome and annotation measurements
Let's get into some computational biology! We will start by obtaining some
basic annotation statistics regarding what is in the database. 

```bash
# return to the main project directory if you have not
cd <PROJECT_DIR>/mycotools_202405

# obtain annotation statistics from the primary MTDB
annotationStats $(mtdb) > annotation_stats.tsv
```

Now you can open that `.tsv` file in whatever spreadsheet software/command-line
text editor you prefer. We can also run a similar analysis for the assemblies:

```bash
assemblyStats $(mtdb) > assembly_stats.tsv
```

Both `annotationStats` and `assemblyStats` can be ran referencing a single
genome `.gff3` or `.fna` respectively IF they have been curated into an MTDB.

<br />

## Circumscribe genes into homology groups and identify single-copy orthologs
It is often useful to group genes into homology groups for downstream analyses.
Single-copy orthologs (SCOs) are a particular type of homology group that is
useful for phylogenetic reconstruction because SCO evolution is assumed to be
congruent with the evolution of the species. A robust automated method of 
SCO determination is OrthoFinder,
but we can more swiftly circumscribe homology groups and identify putative SCOs
using a faster algorithm, `MMseqs cluster`, implemented in `db2hgs` of the
Mycotools suite. Let's run a phylogenetic analysis on a subset of our genomes.

```bash
# extract the lineages we want
mtdb e -l Ustilago >> tree.mtdb
mtdb e -l Puccinia >> tree.mtdb
mtdb e -l Trichoderma >> tree.mtdb
mtdb e -l Psilocybe >> tree.mtdb

db2hgs -d tree.mtdb
```

<br />

## Reconstruct a multigene species phylogeny
A phylogenomic tree is a phylogenetic analysis of a group of genomes that
incorporates an arbitrary minimum number of genes. There are multiple methods
for reconstructing a phylogenomic tree, and Mycotools implements the multigene
partition method. This method essentially determines the best-fit evolutionary
model for each SCO, concatenates an alignment of all SCOs in the dataset, and
then reconstructs a maximum-likelihood phylogenomic tree by applying the best-fit evolutionary
model to each individual gene partion. A general rule-of-thumb is
the more tractable genes we have to construct a phylogenomic tree, the better.
We have quite a few, and it would be ideal to use them all!... but we don't
have time, so let's work with a subset by copying them to a new folder:

```bash
# make a folder for single copy genes
mkdir sco_202405

# copy the top three SCOs
for i in $(ls db2hgs_<YYYYmmdd>/single_copy_genes/ | head -3)
  do cp db2hgs_<YYYYmmdd>/single_copy_genes/$i sco_202405/
done

# run the tree building pipeline
fa2tree -i sco_202405/ --partition
```

When complete, we will open the `concatenated.nex.contree` file in the
resulting `fa2tree_<YYYYmmdd>` directory in FigTree, which is the consensus
tree with 1000 ultrafast bootstrap replicates. 

What you will note is that the tips are labeled with the ome code - but we
probably want to see the actual genus, species, and strain names, right?! Let's
convert the phylogenomic tree from genome code tips to full names:

```bash
ome2name fa2tree_<YYYYmmdd>/concatenated.nex.contree o \
  > fa2tree_<YYYYmmdd>/full_name.newick
```

Go ahead and open this one in FigTree, and let's glance at how well supported
this phylogeny is.

<br />

## Reconstruct the evolution of a gene of interest
Horizontal transfer, gene duplication,
convergence, and gene loss all lead to discordance between the species and
gene evolution. So to study a gene's evolution, we have to reconstruct a
phylogeny of it - let's study the evolution a gene associated with nitrate
assimilation in our dataset.

First, we need to identify homologs of this gene across our database. We will
do this by implementing a BLAST search of the protein sequence. We obtain the
protein sequence using a handy command, `acc2fa`.

```bash
# extract the protein accession of interest
acc2fa -a ustbro1_1795 > ustbro1_1795.faa

# run a blast search on this gene against the primary MTDB
db2search -a blastp -q ustbro1_1795.faa -e 2
```

With the BLAST results in hand, we can now build a phylogeny of the outputted
`.fasta` file of homologs. First, we will move the phylogenomic analysis to a
different directory so `fa2tree` doesn't output to the same place because the
date is the same. This time we will use FastTree to expedite tree
building at the cost of some quality:

```bash
# move the phylogenomic directory
mv fa2tree_<YYYYmmdd> phylogenomic_<YYYYmmdd>/

# run the single gene phylo
fa2tree -i db2search_<YYYYmmdd>/fastas/ustbro1_1795.search.fa -f
```

Now, we can view this tree in FigTree.

<br />

## Visualize the evolution of a gene cluster/syntenic locus
Oftentimes, phenotypes are derived from multiple genes that are proximal to one
another in the linear genome. Therefore, studying the evolution of a locus from 
both a gene-by-gene and locus-wide (shared synteny) perspective can be
important for understanding the evolution a phenotype.

For this example, we will look at the nitrate assimilation gene cluster, which
encodes three genes that import extracellular nitrate and reduce it
to intracellular ammonium as a nitrogen source.

The Cluster Reconstruction and phylogenetic Analysis Pipeline (CRAP)
facilitates swiftly analyzing the evolution of a gene cluster. We can 
implement this pipeline by first extracting a representative nitrate assimilation 
locus, then inputting it into the CRAP pipeline:

```bash
# extract a locus of interest, and store in a file
acc2locus -a ustbro1_1795 -p 1 > nitrate_cluster.txt

# run the CRAP analysis
crap -q nitrate_cluster.txt -s blastp
```

Once complete, we will see `.svg` files in the output directory,
`crap_<YYYYmmdd>` that contain phylogenies of the nitrate assimilation genes,
with synteny diagrams of the individual genes extracted applied to each tip.

<br />

## Create GenBanks of loci for a visually appealing synteny analysis
We can also make more visually appealing synteny diagrams that depict percent
identity through Clinker. I built the CRAP pipeline to output the most similar
loci to the query sequence, and we can invoke this option by having CRAP rerun
off the initial results referencing the initial output directory. This is the
general way for resuming Mycotools scripts.

```bash
# rerun a Mycotools command referencing a previous output directory
crap -q nitrate_cluster.txt -s blastp --loci -o crap_<YYYYmmdd>
```

There should now be a subdirectory, `loci/`, that contains the most similar
loci, and we want to generate GenBanks of these loci to run Clinker. So let's
make a directory in our main project folder:

```bash
mkdir clinker_<YYYYmmdd>
```

Then generate GenBanks of each locus file using some basic BASH scripting:

```bash
for i in crap_<YYYYmmdd>/loci/*txt; do o=$(basename ${i} .txt); acc2gbk -i ${i}
> clinker_<YYYYmmdd>/${o}.gbk; done
```

We can now run Clinker on this dataset:

```bash
cd clinker_<YYYYmmdd>
clinker ./ -p nitrate_assimilation.html
```

<br />

## Integrate external software using Mycotools data extraction
We previously circumscribed homologous gene groups using a script that
implements `MMseqs2 cluster`, which is a fairly crude method. OrthoFinder is a
much more robust method for circumscribing genes into orthologous
gene groups (orthogroups) that are approximately set with respect to the most
recent common ancestral genome of the dataset. OrthoFinder is not built into Mycotools
yet, though we can implement it - and other external scripts - by acquiring the
necessary files.

```bash
# make your OrthoFinder analysis directory
mkdir orthofinder_<YYYYmmdd>/
cd orthofinder_<YYYYmmdd>/

# extract the proteomes you need for an analysis, referencing a MycotoolsDB we
# extracted previously
db2files -p -d ../tree.mtdb
```

This will generate a folder of proteome fastas, `faa/`, that we can now
reference. Note this folder contains "symlinked" files, which are essentially
empty links to the original file. Therefore, editing these files will edit the
initial files in the database. If you want to make hard copies of the files,
append the `--hard` flag to the `db2files` command.

```bash
orthofinder -f faa/
```

If your dataset is small enough then I recommend implementing
OrthoFinder over `db2hgs`

<br /><br /><br />

## REFERENCE SOFTWARE

#### BioPython:

Chapman, B. & Chang, J. Biopython: Python Tools for Computational Biology. SIGBIO Newsl. 20, 15–19 (2000).

#### BLAST:

Ye, J., McGinnis, S. & Madden, T. L. BLAST: improvements for better sequence analysis. Nucleic Acids Research 34, W6–W9 (2006).

#### Clinker:

clinker & clustermap.js: automatic generation of gene cluster comparison figures | Bioinformatics | Oxford Academic. https://academic.oup.com/bioinformatics/article/37/16/2473/6103786.

#### ClipKIT:

Steenwyk, J. L., Iii, T. J. B., Li, Y., Shen, X.-X. & Rokas, A. ClipKIT: A multiple sequence alignment trimming software for accurate phylogenomic inference. PLOS Biology 18, e3001007 (2020).

#### CRAP:

Slot, J. C. & Rokas, A. Horizontal Transfer of a Large and Highly Toxic Secondary Metabolic Gene Cluster between Fungi. Current Biology 21, 134–139 (2011).

#### FastTree:

Price, M. N., Dehal, P. S. & Arkin, A. P. FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. PLoS ONE 5, e9490 (2010).

#### IQ-TREE:

Minh, B. Q. et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. Molecular Biology and Evolution 37, 1530–1534 (2020).

#### MAFFT:

Katoh, K., Misawa, K., Kuma, K. & Miyata, T. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res 30, 3059–3066 (2002).

#### MMseqs2:

Steinegger, M. & Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35, 1026–1028 (2017).

#### OrthoFinder:

Emms, D. M. & Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biology 20, 238 (2019).

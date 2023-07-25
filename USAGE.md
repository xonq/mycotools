![Mycotools](https://gitlab.com/xonq/mycotools/-/raw/master/misc/pictogo.png)

# OVERVIEW
Mycotools is a comparative genomics software suite centered around a curated
database (MycotoolsDB/MTDB) of genomic data. This guide references Mycotools scripts that enable
high throughput pipelining and routine genomic analysis. Mycotools is *HEAVILY*
dependent on the MycotoolsDB, so there is no guarantee that these scripts will
work with external files. 

Why the dependency on MycotoolsDB? Because bioinformatics file formats need systematic, 
uniform curation. [See here for more](https://github.com/xonq/mycotools/blob/master/MTDB.md)

<br /><br />

# USAGE GUIDE
**Table of Contents**

<br />


- **MYCOTOOLSDB TOOLS**
	- [Initializing *de novo* MycotoolsDB](https://github.com/xonq/mycotools/blob/master/USAGE.md#de-novo-initialization)
        - [Initializing reference MycotoolsDB](https://github.com/xonq/mycotools/blob/master/USAGE.md#reference-initialization)
	- [Updating MycotoolsDB](https://github.com/xonq/mycotools/blob/master/USAGE.md#updating)
	- [Connecting to the database](https://github.com/xonq/mycotools/blob/master/USAGE.md#interfacing)
	- [Managing the database](https://github.com/xonq/mycotools/blob/master/USAGE.md#managing)
	- [Querying the database](https://github.com/xonq/mycotools/blob/master/USAGE.md#querying)
	- [Creating modular databases](https://github.com/xonq/mycotools/blob/master/USAGE.md#creating-modular-databases)
	- [Acquiring database files / file paths](https://github.com/xonq/mycotools/blob/master/USAGE.md#acquiring-database-files)
	- [Adding local genomes to the database](https://github.com/xonq/mycotools/blob/master/USAGE.md#adding-local-genomes)
	- [Substitute organism name for MycotoolsDB organism code](https://github.com/xonq/mycotools/blob/master/USAGE.md#ome2namepy)

<br />


- **SEQUENCE DATA**
	- [Downloading from NCBI / JGI](https://github.com/xonq/mycotools/blob/master/USAGE.md#downloading-files)
	- [Sequence data statistics](https://github.com/xonq/mycotools/blob/master/USAGE.md#sequence-data-statistics)
	- [Grabbing accessions](https://github.com/xonq/mycotools/blob/master/USAGE.md#grab-accessions)
	- [Grabbing full GenBank](https://github.com/xonq/mycotools/blob/master/USAGE.md#full-genome-genbanks)
	- [Extract fasta coordinates](https://github.com/xonq/mycotools/blob/master/USAGE.md#fasta-coordinates)
	- [Grabbing loci](https://github.com/xonq/mycotools/blob/master/USAGE.md#grab-loci)
	- [GFF to sequence](https://github.com/xonq/mycotools/blob/master/USAGE.md#gene-coordinates-to-sequences)
	- [Curating annotation](https://github.com/xonq/mycotools/blob/master/USAGE.md#curate-annotation)
	- [Adding corrected gene models](https://github.com/xonq/mycotools/blob/master/USAGE.md#adding-corrected-gene-models)
	- [Visualizing loci](https://github.com/xonq/mycotools/blob/master/USAGE.md#visualizing-loci)
<br />


- **EVOLUTIONARY ANALYSES**
	- [MycotoolsDB BLAST/HMM](https://github.com/xonq/mycotools/blob/master/USAGE.md#homolog-search-mycotoolsdb)
	- [Fasta to tree](https://github.com/xonq/mycotools/blob/master/USAGE.md#tree-building)
	- [Sequence clustering](https://github.com/xonq/mycotools/blob/master/USAGE.md#sequence-clustering)
	- [Gene cluster phylogenetic analysis](https://github.com/xonq/mycotools/blob/master/USAGE.md#crappy)


<br />

- **MYCOTOOLS PIPELINES**
    - [Mycotools pipelining principles](https://github.com/xonq/mycotools/blob/master/USAGE.md#shell-pipelining-with-mycotools)
	- [Phylogenetic analysis](https://github.com/xonq/mycotools/blob/master/USAGE.md#phylogenetic-analysis)


<br />

![Analysis_examples](https://github.com/xonq/mycotools/blob/master/misc/examples.png?raw=true)

---



<br /><br /><br />


# MYCOTOOLSDB
The following scripts interface with and manipulate MycotoolsDB (MTDB) `.mtdb`
files. To learn more about MycotoolsDB and the `.mtdb` format standard, refer to
[this guide](https://github.com/xonq/mycotools/blob/master/MTDB.md).

<br /><br />

## Initialization
### mtdb update
`mtdb update` is for initializing and building the primary MTDB; to interface with
an established MTDB, see
[interfacing](https://github.com/xonq/mycotools/blob/master/USAGE.md#interfacing)

### De novo initialization
#### FUNGI
To initialize a curated database of all NCBI and MycoCosm (JGI) fungal genomes; note, 
JGI downloading is limited to one file per minute, so initialization will take multiple days
and can be resumed with `-r YYYYmmdd`:

```bash
mtdb update -i <INIT DIRECTORY>
```

#### PROKARYOTES
To initialize a curated database of all NCBI prokaryotic genomes (please note
this is in alpha-testing and needs to improve vectorization to scale more efficiently; additionally, there are 100,000s prokaryote genomes to download and this will take several days-weeks):
```bash
mtdb update -i <INIT_DIRECTORY> -p
```

<br />

### Reference initialization
A MycotoolsDB can be initialized referencing an external `.mtdb` file, e.g. to reproduce an analysis using another dataset. Please note this will only assimilate NCBI and JGI genomes from the reference `.mtdb`.

If you are currently linked to an existing MycotoolsDB, unlink via:
```bash
mtdb -u
```

To initialize a MycotoolsDB from a reference, appending necessary arguments:
```bash
mtdb update -i <INIT_DIR> -r <REF.mtdb>
```

If successful, a new MycotoolsDB will be initialized in `<INIT_DIR>`; to link back to any previously established MycotoolsDBs see how to [interface](https://github.com/xonq/mycotools/blob/master/USAGE.md#interfacing). 

<br /><br />

## Updating
To update MycotoolsDB:

```bash
mtdb update -u
```

<br /><br />

## Interfacing
### mtdb
`mtdb` is the MycotoolsDB central utility. It initializes interfacing with an established primary
database or just prints the path of the primary database. MycotoolsDBs are labelled `YYYYmmdd.mtdb`.
```bash
mtdb
/home/xonq/mtdb/mtdb/20210125.mtdb
```

To add interfacing with a fungal/prokaryote primary MTDB:
```bash
mtdb -i <DATABASE_BASE_DIRECTORY>
```

To switch between established interfaces:
```bash
mtdb -f
mtdb -p
```

NOTE: only one MycotoolsDB of each type will be stored. Multiple databases of the same type require reiniterfacing via `mtdb -i <PATH>`.

<br /><br />

## Managing
### mtdb manage
To encrypt a local copy of your NCBI and JGI passwords for fast access, run:
```bash
mtdb manage -p
```

This is necessary to submit jobs for scripts that require passwords. Such scripts 
(like `mtdb update`) will need to receive a password from stdin e.g.:

```bash
read -s PASSWORD
<type password>
printf '{"username":"myname","password":"%s"}' $PASSWORD | mtdb update -u
```

<br /><br />

## Querying
`mtdb` can query ome codes. To obtain the whole row for an ome code, simply:
```bash
mtdb <OME>
```

To obtain the specific file path for an ome file, append `.fna` (assembly),
`.faa` (proteome), or `.gff3` (gene coordinates) to the ome:
```bash
mtdb <OME>.<EXTENSION>
```

Because you can obtain PATHs using `mtdb`, you can use basic bash functionality to work with the output, 
e.g. to open in a text editor or to grep the file:
```bash
vim $(mtdb)
grep 'Psilocybe' $(mtdb)
```

NOTE: `grep` may yield non-specific taxonomy results, e.g. `grep Athelia
$(mtdb)` will not only yield the genus *Athelia*, but also all members of the
family Atheliaceae. See
[mtdb extract](https://github.com/xonq/mycotools/blob/master/USAGE.md#creating-modular-databases)
to extract based on taxonomy.

<br /><br />

## Creating modular databases
### mtdb extract
`mtdb extract` subsets lineages of interest from the primary MTDB. By default, this script will ignore use-restricted data - add `-n` to include use-restricted data if you are aware of and will respect the limitations of use-restricted data. Run `mtdb extract -h` to see all options.

e.g. grab a database of a taxonomic order: 
```bash
mtdb extract -l Atheliales > atheliales.mtdb
```

grab all NCBI Aspergilli accessions: 
```bash
mtdb extract -s ncbi -l aspergillus > aspergillus.mtdb_ncbi
``` 

grab a list of orders from a file:
```bash
mtdb extract -ll <TAX_FILE> > taxa.mtdb
```

<br /><br />


## Acquiring database files
### db2files
Inputs a `.mtdb` file (by default uses the primary database), then symlinks the selected file types, hard copies the files, or prints their PATHs. A symlink is simply creating a placeholder file that links to the database file... this way it does not take up additional storage space like a hard copy does. However, editing symlinks will edit the original file, so *only hard copy `--hard` if you need to edit the files*.

Let's say you want protein data from organisms in one family. First, you should extract a database of organisms you want:
```bash
mtdb extract -l Atheliaceae > atheliaceae.mtdb
```

Then, run `db2files` to copy the protein fastas into the current directory (call `-h` to see all options):
```bash
db2files -d atheliaceae.mtdb -p 
```

Alternatively, if you just need the paths (links) to these files, simply run:
```bash
db2files -d atheliaceae.mtdb -p --print
```

<br /><br />

## Adding local genomes
### mtdb predb2mtdb
To add in-house annotations `mtdb predb2mtdb` will input your genome and
metadata, curate, and prepare a database file to add to the database. The
administrator will then take your database and add it to the primary MTDB.

First, generate a predb spreadsheet:
```bash
mtdb predb2mtdb > predb.tsv
```

The resulting `predb.tsv` can be filled in via spreadsheet software and
exported as a tab delimited `.tsv`. Alternatively, use a plain text editor and
separate by tabs. De novo annotations produced by Funannotate/Orthofiller must
be filled in as "new" for the genomeSource column; *annotations directly
derived from NCBI/JGI data need to be specified in genomeSource. Updates to
existing database entries may not be automatically detected using the metadata,
so please explicitly enter the current MTDB ome code you are updating in the
`previous_ome` column.

Next, generate a MycotoolsDB file from your completed predb, and notify your
database administrator that it is ready for integration:
```bash
mtdb predb2mtdb <PREDB.TSV>
```

Finally, inform your database administrator that the database is ready for
submission to the primary MTDB. Administrators will execute:
```bash
mtdb update -a <PREDB2DB.mtdb>
```

<br /><br />


## Other MycotoolsDB scripts
### ome2name
Substitutes MycotoolsDB organism code names (e.g. `fusgra1`) for taxonomic information (e.g. Fusarium_graminearum_var._XYZ).

e.g. to substitute ome for genus species and strain: `ome2name <INPUT> oa`
```bash
ome2name -h
USAGE: ome2name <INPUTFILE> | ome2name <INPUTFILE> [MYCOTOOLSDB] asvg*&
DEFAULTS: primary db, see script for default forbidden characters
Input file to regex sub omes with their name.
optional MycotoolsDB, string of forbidden characters
"o" no ome | "g" no genus | "s" no species | "v" no strain | "a" no alternative ome
```

<br /><br /><br />

# SEQUENCE DATA TOOLS
## Sequence data statistics
### assemblyStats / annotationStats
```bash
assemblyStats <ASSEMBLY.fa>
annotationStats <ANNOTATION.gff3>
```

To obtain a table of annotation statistics, [create a mycotoolsDB](https://github.com/xonq/mycotools/blob/master/USAGE.md#creating-modular-databases) file with the organisms of interest and run:
```bash
assemblyStats <MYCOTOOLSDB.mtdb>
annotationStats <MYCOTOOLSDB.mtdb>
```

If you want to route the output to a file, simply redirect output by appending ` > <OUTPUTFILE>` to the command, or add an output file as the second argument


<br /><br />


## Downloading files
### jgiDwnld / ncbiDwnld
These scripts input a MycotoolsDB or can be manually made as shown at the bottom of this section. 

Say you want to grab transcript information from a genus, *Aspergillus*. First, extract entries in the database that are within *Aspergillus*:
```bash
mtdb extract -l aspergillus > aspergillus.mtdb_ncbi
```

If there are organisms you do not want in the extracted `.mtdb`s, delete their line(s) in the file. Next call `jgiDwnld -h` or `ncbiDwnld -h` to find the flags necessary to download the files you want. To download transcript data (and EST data for JGI) in your current directory:
```bash
jgiDwnld -i aspergillus.mtdb_jgi -t -e
ncbiDwnld -i aspergillus.mtdb_ncbi -t
```

To unzip all the files, run `gunzip <FILETYPE>/*.gz`. 
To submit as a job (not recommended), you must create an encrypted MycotoolsDB passkey using [`mtdb manage`](https://github.com/xonq/mycotools/blob/master/USAGE.md#managing) and pass the password to stdin to these scripts.

<br />

These scripts can input assembly accessions (NCBI) or genome codes (JGI). The column must have the appropriate header ('assembly_acc' or 'assembly_acc'):

`jgiGenomeCodes.txt`
```
assembly_acc
Abobi1
Absrep1
Acain1
```

`ncbiBioSamples.txt`
```
assembly_acc
SAMN15352002
SAMN06205226
SAMN02744098
```

<br />

You can download NCBI SRA's after downloading NCBI's SRA tools and making sure `fastq-dump` is in your PATH:

```bash
ncbiDwnld --sra -i <REFERENCE>
```

You can create a file with SRA ID's or BioProject, etc. Basically any query that is unique and sufficient to acquire the SRRs of interest. For paired-end reads, include `-pe`.

<br /><br />


## Grab accessions
### acc2fa / acc2gff / acc2gbk
All Mycotools accessions - assembly or protein - are `<ome>_<acc>` where "ome"
is the MTDB genome codename and "acc" is the accession.
Let's say you want to query *Panaeolus cyanescens'* (ome: pancya1) PsiD and the NCBI accession is "PPQ80975.1"

Prepend the ome codename to "_" and the accession and grab it from the MTDB:
```bash
acc2gff -a pancya1_PPQ80975.1
acc2fa -a pancya1_PPQ80975.1
acc2gbk -a pancya1_PPQ80975.1
```

To output to file, simply append `> <FILENAME>` to your command. 

If you have a list of accessions, create an input file with the accessions separated by new lines then run:
```bash
acc2gff -i <INPUTFILE>
acc2fa -i <INPUTFILE>
acc2gbk -i <INPUTFILE>
```

Alternatively, input accessions from stdin:
```bash
echo "pancya1_PPQ80975.1" | acc2gbk -a -
```

<br />

### Full genome GenBanks
Full genome GenBanks (GBKs) are not available by default to save space.
To retrieve a full genome GBK, simply append `--full` to the `acc2gbk`
command and change `-a` from an accession to an ome code.

```bash
acc2gbk -f -a athter2.2 > athter2.2.full.gbk
```

<br />

### Grab loci
### acc2locus
Grab loci the same as above within a set number of genes or nucleotides (`-n`) plus or minus:

list proximal +/- 5Kb from accession to standard out; NOTE `-n` for nucleotide:
```bash
acc2locus -a fibpsy1_906341 -p 5000 -n
fibpsy1_809145
fibpsy1_923701
fibpsy1_771516
fibpsy1_906341
fibpsy1_719531
fibpsy1_846242
fibpsy1_138
```

list proximal +/- 5 genes to standard out; NOTE no `-n`:
```bash
acc2locus -a fibpsy1_906341 -p 5
fibpsy1_880711
fibpsy1_846234
fibpsy1_809145
fibpsy1_923701
fibpsy1_771516
fibpsy1_906341
fibpsy1_719531
fibpsy1_846242
fibpsy1_138
fibpsy1_942299
fibpsy1_906343
```

grab genes between two accessions; NOTE `-b` for between:
```bash
acc2locus -a "psicub1_30114 psicub1_87205" -b
psicub1_30114
psicub1_72370
psicub1_30121
psicub1_30049
psicub1_72373
psicub1_87205
```

Generate a genbank from the locus (useful for `clinker`):
```bash
acc2locus -a <OME>_<ACC> -p 5 | acc2gbk -a -
```

Output GFFs and protein fastas for the locus:
```bash
acc2locus -a <OME>_<ACC> -p 5 --ome
```

<br /><br />

## Fasta coordinates
### coords2fa

```bash
coords2fa -h

Input nucleotide fasta/tsv input, extract coordinates
coords2fa <FA> <SEQID> <START_COORD> <END_COORD> <STRAND_SENSE>
Extract full sequence from sense strand: coords2fa test.fna scaffold_20 0 -1
Extract coordinates from antisense strand: coords2fa test.fna scaffold_20 69 420 -


Bulk extraction tab delimitted row format:
#fasta_path   sequence_id     start_coordinate        end_coordinate  strand_sense

coords2fa coords.tsv
```

<br /><br />

## Gene coordinates to sequences
### gff2seq
`gff2seq` will extract the nucleotide or amino acid sequences associated with a gene coordinates `gff` file. Accepts optional flanking nucleotide plus/minus to extract from flanks of coordinates. `gff2seq` will extract flanks independently for each sequence ID (column 1) within the `gff`. Coding or noncoding regions can be specified. If `--intergenic` is called then only the first and last gene for each fasta
sequence are considered.

e.g. extract nucleotide sequences and 1 kilobase flanks and noncoding regions within the following genes:

```bash
gff2seq -g <.GFF3> -a <.FNA> -nc -pm 1000 -n
```


<br /><br />

## Curate annotation
Please note that MTDB curates annotations and assemblies prior to adding to the database via
`mtdb predb2mtdb`. If your intention is to add your data to the database, please see `mtdb predb2mtdb`.
If your data can be curated, it will be submitted to one of the following scripts during 
`mtdb predb2mtdb`. These scripts are not in your PATH by default and are useful for testing 
compatibility with MTDB. Only complete genome annotations are acceptable. 
MTDB strives to incorporate the most common annotation software -
please raise an issue if your complete genome annotation - from a common software - is not compatible.

### gtf2gff3
This script will convert OrthoFiller `.gtf` to `.gff3`, rename headers sequentially, and 
optionally adds Alias

```bash
gtf2gff3 -g <ORTHOFILLER>/results/results.gtf -f <ASSEMBLY> -p <OME>
```

### gff2gff3
Designed for JGI gff2 legacy formatted gffs
```bash
gff2gff3 -o <OME> -j <JGI_ome> -i <GFF2>
```

### curGFF3 
`curGFF3` is tested with, Funannotate, JGI, and GenBank `gff3` files
```bash
curGFF3 <GFF3> <OME>
```

<br /><br />

## Adding corrected gene models
### add2gff
This script will add new and corrected gene models to a full genome gff.
Overlapping coordinates can be removed and an update file for the primary MTDB  is optionally generated.

To add corrected genes from an `exonerate`-derived gff to a MycotoolsDB `gff` and prepare a database update:

`add2gff -i <EXONERATE_GFF> -a $(mtdb <OME>.gff3) -u`

NOTE: database administrators will need to finalize the updates to propagate
the data to the primary database via:

```bash
mtdb update -a <ADD2GFF.mtdb>
```

<br /><br />

## Visualizing Loci
### gff2svg
This script will input a .gff3, or new line-delimited list of .gff3 paths, 
and for each contig output an SVG of the locus annotated by function retrieved 
from the 'product=' field. 

<br />

e.g. make an SVG from a GFF3: `gff2svg -g <MY.gff3>`

make SVGs for all GFF3s in a new line delimited list with width set to 20:

```bash
gff2svg -i <LISTOFGFF3.nsv> -o <OUTPUT_DIR> -w 20
```


<br /><br /><br />

# EVOLUTIONARY ANALYSIS TOOLS
## Homolog search MycotoolsDB
### db2search
`db2search` will execute `hmmer`, `blast`, `diamond`, or `mmseqs`, query an input fasta, and output a results fasta for each accession in the query.

<br /><br />


## Tree Building
### fa2tree
`fa2tree` will input a fasta file or directory of fasta files, trim, and
generate trees. 

Information on a complete [phylogenetic pipeline](https://github.com/xonq/mycotools/blob/master/USAGE.md#mycotools-pipelines) is elaborated below.

<br />

e.g. Swiftly construct a fastree from a fasta
```bash
fa2tree -i <FASTA>.fa --fast
```

Prepare a robust IQ-TREE for slurm job submission
```bash
fa2tree -i <FASTA> -s -A PAS1046
```

Construct a multigene phylogeny with independent evolutionary models for each gene
```bash
fa2tree -i <FASTA_DIR> -p
```

The final IQ-TREE file you want is `*.contree`; fastree is a `.treefile`
View your trees using [FigTree](https://github.com/rambaut/figtree/releases).

<br /><br />


## Sequence clustering
### fa2clus
Some gene families (e.g. P450s) have many highly identical homologs, which
is problematic for conducting phylogenetic analysis and manipulating these large datasets. 
`fa2clus` invokes sequence clustering algorithms to systematically truncate your dataset
without constructing a phylogeny. `fa2clus` optionally implements an automated iterative 
approach to obtaining a cluster of minimum - maximum size with the gene of interest.

For hierarchical agglomerative clustering: `fa2clus` will either take a `fasta` and generate a distance matrix using 
`usearch calc_distmx` by default or the % identity of `diamond` alignments.
Then, cluster sequences via hierarchical agglomerative clustering and output a
`.clus` file of cluster assignments and `.newick` dendrogram. 

Please note that sequence similarity clustering is using a heuristic, sequence
similarity, as a proxy of relatedness. Therefore, sequence similarity
clustering may remove closely related sequences in cases where the inputted number of
genes is much greater than the maximum cluster constraint and/or if there are
not highly similar sequences around the gene of interest

To run `mmseqs cluster` on a faster with minimum 20% query coverage and minimum 30% AA identity

```bash
fa2clus -f <FASTA>.fa -m 0.2 -x 0.3
```

Iteratively cluster until a cluster size of 50-200 genes is achieved:
```bash
fa2clus -f <FASTA> -m 0.2 -x 0.3 --iterative <FOCAL_GENE> --minseq 50 --maxseq 200
```

<br /><br />


## Gene Cluster Reconstrunction and Phylogenetic Analysis (CRAP)
### crap

<img align="left"
src="https://github.com/xonq/mycotools/blob/master/misc/crap_example.png?raw=true"
alt="Extracted clade of CRAP pipeline" height="450" width="578">

CRAP, adopted and expanded from [Slot & Rokas implementation](https://doi.org/10.1016/j.cub.2010.12.020),
reconstructs and visualizes gene cluster phylogenies to study gene cluster
evolution on a gene-by-gene basis. CRAP will: 

1) input a cluster query and use a search algorithm (BLAST/mmseqs/Diamond)
or orthogroup-based approach to find homologs in the MycotoolsDB

2) implement sequence similarity clustering to truncate the sequence set 
and detect outgroups. NOTE: sequence similarity clustering can sometimes remove close homologs
if there are too many sequences analyzed relative to the set maximum sequences. 

3) construct phylogenies of each query sequence

4) map locus synteny diagrams onto the tips of the phylogenies.

`crap` can operate on a query of MycotoolsDB accessions or a standalone
multifasta input of external accessions. Following homolog acquisition,
`crap` will submit each set of hits for tree building or sequence similarity
clustering if the number of sequences exceeds the inputted maximum 
(`-m`).

By default, `crap` will construct trees using fasttree.
Alternatively, CRAP can construct a robust IQ-TREE with 1000 boostrap iteration support values
by specifying `-i`/`--iqtree`.

<br />

To search an extracted sub-MycotoolsDB using `blastp` and create phylogenies with `fasttree`:
```bash
mtdb extract --lineage Basidiomycota > basi.mtdb
crap -q <QUERYGENES> -d basi.mtdb -s blastp --bitscore 40 --cpu 12
```

<br /><br /><br />

# Mycotools Pipelines
## Shell Pipelining with Mycotools
Mycotools is designed to enable pipelining in Linux shells as well as python.

For scripts such as `acc2fa/acc2gbk/acc2gff/acc2locus`, the input is an accession or set of
accessions, and can be piped in via standard input. *All scripts that accept
standard input (stdin) will require "-" as the input argument.* 
Stdin access allows you can to chain commands and swiftly generate the necessary
input for downstream analysis.

Say you want to see if an accession of interest is part of a gene cluster. Run a CRAP
around accessions within 20 kb +/- your gene of interest:
```
acc2locus -a athter2_3 -n -p 20000 | crap -q - <ARGS>
```

Create a file of query genes from the CRAP output, extract locus accessions for each gene,
and output a genbank for each locus:

```bash
for acc in $(cat <INPUT_FILE>)
do 
  acc2locus -a $acc -p 20000 -n | acc2gbk -a - > $acc.locus.gbk 
done
```

<br />

## Phylogenetic analysis

Despite the benefits of increased sampling, there are two prominent problems
large samples create in phylogenetic analysis: 1) large alignments
lose resolution because they must consider many sequences and 2) analysis
time increases with sample size. 
For a small set of genes, such as ITS, one can usually assume that the gene family is 
conserved, so the dataset can be cut down to closely 
related organisms. For most other genes, it is not valid to assume the gene family
is conserved because horizontal transfer can lead to unexpected distributions.

It is thus sometimes necessary to balance alignment resolution and computational
tractability by systematically truncating the dataset into subset of gene homologs. 
This is accomplished in Mycotools by iteratively constructing phylogenies, 
identifying and extracting clades of homologs,
and repeating until a manageable tree is obtained. On its own, this analysis requires 
elaborate integration of multiple independent softwares, but Mycotools takes care of most of this.

<br />

### Example 1:

Acquire a set of homologs for a gene family of interest by BLASTing a representative
query protein sequence. 

1. extract a database of published sequences, or use other arguments to extract
other organisms of interest

```bash
mtdb extract > pub.mtdb
```

<br />

2. obtain gene homologs using `db2search` with an e-value threshold of 10<sup>-2</sup>:

```bash
db2search -d pub.mtdb -e 2 -q <QUERY>
```

3. If there are more than 1000 homologs, truncate the resulting homologs around a gene of 
interest between 50 and 500 genes:

```bash
fa2clus -f <FASTA>.fa --min_seq 50 --max_seq 500 -i <FOCAL_GENE>
```

4. If there are less than 1000 homologs, construct a fastree, `-f`, if there are many samples,
or remove that parameter for a robust IQ-TREE.

```bash
fa2tree -i <FASTA>.fa 
```

5. If the phylogeny can be improved, extract a highly supported node from 4 and reconstruct,
or systematically truncate the dataset with 3 and reconstruct.

<img align="right" src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/ablogo.png">

<br /><br /><br /><br /><br /><br /><br /><br /><br />

## TODO
- [ ] add ncbiAcc2fa
- [ ] separate better/organize into wiki
- [x] document clinker pipeline
- [x] add images of pipeline output
- [ ] fa2hmmer2fa annotation
- [x] update crap information
- [ ] eggnog to synteny diagram
- [ ] overview of mycotools scripts functions, e.g. resume on -o, output standards
- [x] update phylogenetic pipeline with fa2clus renovations
- [x] add partition analysis tutorial
- [ ] discuss output in crap
- [ ] add bioreform
- [ ] convert to wiki
- [ ] add caveats to outputs of fa2clus, crap, etc

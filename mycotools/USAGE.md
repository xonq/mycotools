# INDEX

<br />


- **MYCOTOOLSDB TOOLS**
	- [Initializing MycotoolsDB](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#initialization)
	- [Updating MycotoolsDB](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#updating)
	- [Interfacing with the database](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#mycotools-db)
	- [Creating modular databases](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#creating-modular-databases)
	- [Acquiring database files / file paths](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#acquiring-database-files)
	- [Adding local genomes to the database](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#adding-local-genomes)
	- [Substitute organism name for MycotoolsDB organism code](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#ome2name.py)

<br />


- **SEQUENCE DATA**
	- [Downloading from NCBI / JGI](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#downloading-files)
	- [Sequence data statistics](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#sequence-data-statistics)
	- [Grabbing accessions](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#grab-accessions)
	- [Extract fasta coordinates](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#fasta-coordinates)
	- [GFF to sequence](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#gene-coordinates-to-sequences)
	- [Grabbing loci](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#grab-loci)
	- [Visualizing loci](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#visualizing-loci)
	- [Curating annotation](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#curate-annotation)

<br />


- **EVOLUTIONARY ANALYSES**
	- [MycotoolsDB BLAST](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#blast-mycotoolsdb)
	- [MycotoolsDB hidden markov model search](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#hmmsearch-mycoDB)
	- [Fasta to tree](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#tree-building)
	- [Hierarchical agglomerative clustering](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#hierarchical-agglomerative-clustering)

<br />

- **MYCOTOOLS PIPELINES**
	- [Phylogenetic analysis](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#phylogenetic-analysis)
	- [Cluster Reconstruction and Phylogenetic Analysis (CRAP)](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#crap.py)


---

<br /><br /><br />


# MYCOTOOLS DB
The following are information regarding scripts that manipulate Mycotools .mtdb files. To learn more about Mycotools .mtdb, you may refer to the master [MycotoolsDB repository](https://gitlab.com/xonq/mycotoolsdb/-/README.md) (currently unavailable).

<br /><br />

## Initialization
### updateDB.py
To initialize a curated database of all NCBI and MycoCosm (JGI) genomes, run the following

```bash
updateDB.py -i <INIT DIRECTORY>
```

<br /><br />

## Updating
To update MycotoolsDB:

```bash
updateDB.py -u
```

<br /><br />

## Interfacing with the master database
### mycodb
`mycodb` is a utility that integrates with the master database or just prints the path of the master database, which can then be used with other shell commands. MycotoolsDBs are labelled `YYYYmmdd.mtdb`.
```bash
mycodb
/home/xonq/mycodb/mycodb/20210125.mtdb
```

If you want to use the path of the master database you can use basic bash functionality to work with the output, 
e.g. to open in a text editor or to grep the file:
```bash
vim $(mycodb)
grep 'Psilocybe' $(mycodb)
```

<br /><br />


## Creating modular databases
### extractDB.py
If you are only interested in a subset of lineages in the master mycotoolsDB, then extract the portion you want, run `extractDB.py -h` to see all options:

e.g. grab a database of a taxonomic order: 
```bash
extractDB.py -l Atheliales -r order > atheliales.mtdb
```

grab all NCBI Aspergilli accessions: 
```bash
extractDB.py -s ncbi -l aspergillus -r genus > aspergillus.mtdb_ncbi
``` 

grab the inverse of your arguments: 
```bash
extractDB.py -s ncbi -l aspergillus -r genus -i > notAspergullis.mtdb_notNcbi
```

grab a list of orders from a file:
```bash
extractDB.py -ll <TAX_FILE> -r order > taxa.mtdb
```

grab a list of `ome`s in a new line delimited file: 
```bash
extractDB.py ---ome <OME_FILE>
```

```bash
extractDB.py --help
usage: extractDB.py [-h] [-d DATABASE] [-l LINEAGE] [-r RANK] [-s SOURCE] [-n]
                    [-u] [--unique_strains] [-i] [--headers] [-o OUTPUT] [-]
                    [-ol OME] [-ll LINEAGES]

Extracts a MycotoolsDB from arguments. E.g. `extractDB.py -l Atheliaceae -r
family`

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        DEFAULT: masterdb
  -l LINEAGE, --lineage LINEAGE
  -r RANK, --rank RANK  Taxonomy rank
  -s SOURCE, --source SOURCE
                        Data source
  -n, --nonpublished    Include restricted-use
  -u, --unique          Unique species
  --unique_strains
  -i, --inverse         Inverse arguments
  --headers             Include header
  -o OUTPUT, --output OUTPUT
  -, --stdin            Pipe MycotoolsDB from stdin
  -ol OME, --ome OME    File w/list of omes
  -ll LINEAGES, --lineages LINEAGES
                        File w/list of lineages (same rank)
```

<br /><br />


## Acquiring database files
### db2files.py
Inputs a MycotoolsDB `.mtdb` file (by default uses the master database), then creates symlinks of the selected file types, hard copies the files, or prints their PATHs. A symlink is simply creating a placeholder file that links to the database file... this way it does not take up additional storage space like a hard copy does. However, editing symlinks will edit the original file, so *only hard copy `--hard` if you need to edit the files*.

Let's say you want protein data from organisms in one family. First, you should extract a database of organisms you want:
```bash
mkdir pullFiles && cd pullFiles
extractDB.py -r family -l Atheliaceae > atheliaceae.mtdb
```

Then, run `db2files.py` to copy the protein fastas into the current directory (call `-h` to see all options):
```bash
db2files.py -d atheliaceae.mtdb -p 
```

Alternatively, if you just need the paths (links) to these files, simply run:
```bash
db2files.py -d atheliaceae.mtdb -p --print
```

<br /><br />

## Adding local genomes
### predb2db.py
To add in-house annotations predb2db.py will input your genome and metadata, curate, and prepare a database file to add to the database. The manager of the database (Zach for Ohio Supercomputer) will then take your database and add it to the master.

First, generate a predb spreadsheet:
```bash
predb2db.py > predb.tsv
```

The resulting `predb.tsv` can be filled in via spreadsheet software and exported as a tab delimited `.tsv`. Alternatively, use a plain text editor and separate by tabs. De novo annotations produced by Funannotate/Orthofiller must be filled in as "new" for the genomeSource column. If using a plain-text editor, export with UTF-8 formatting and account for blank entries with a tab.

Finally, generate a mycotoolsDB file from your predb, and notify your database manager that it is ready for integration:
```bash
predb2db.py <PREDB.TSV>
```

<br /><br />


## Other MycotoolsDB scripts
### ome2name.py
Substitutes MycotoolsDB organism code names (e.g. `fusgra1`) for taxonomic information (e.g. Fusarium_graminearum_var._XYZ).

e.g. to substitute ome for genus species and strain: `ome2name.py <INPUT> oa`
```bash
ome2name.py -h
USAGE: ome2name.py <INPUTFILE> | ome2name.py <INPUTFILE> [MYCOTOOLSDB] asvg*&
DEFAULTS: master db, see script for default forbidden characters
Input file to regex sub omes with their name.
optional MycotoolsDB, string of forbidden characters
"o" no ome | "g" no genus | "s" no species | "v" no strain | "a" no alternative ome
```

<br /><br /><br />

# SEQUENCE DATA TOOLS
## Sequence data statistics
### assemblyStats.py / annotationStats.py
```bash
assemblyStats.py <ASSEMBLY.fa>
annotationStats.py <ANNOTATION.gff3>
```

To obtain a table of organisms' annotation statistics, [create a mycotoolsDB](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#creating-modular-databases) file with the organisms of interest and run:
```bash
assemblyStats.py <MYCOTOOLSDB.mtdb>
annotationStats.py <MYCOTOOLSDB.mtdb>
```

If you want to route the output to a file, simply redirect output by appending ` > <OUTPUTFILE>` to the command, or add an output file as the second argument


<br /><br />


## Downloading files
### jgiDwnld.py / ncbiDwnld.py
These scripts input a MycotoolsDB or can be manually made as shown at the bottom of this section. 

Say you want to grab a few organisms' transcript information from your genus, *Aspergillus*. First, extract entries in the database that are within *Aspergillus*:
```bash
mkdir dwnldFiles && cd dwnldFiles
extractDB.py -c genus -t aspergillus > aspergillus.mtdb_ncbi
```

If there are organisms you don't want in the extracted `.mtdb`s, just delete their line(s) in the file. Next call `jgiDwnld.py -h` or `ncbiDwnld.py -h` to find the flags necessary to download the files you want. To download transcript data (and EST data for JGI) in your current directory:
```bash
jgiDwnld.py -i aspergillus.mtdb_jgi -t -e
ncbiDwnld.py -i aspergillus.mtdb_ncbi -t
```

These scripts populate with compressed files. To unzip all the files, run `gunzip <FILETYPE>/*.gz`. You will also see log files for the download process. 
To submit as a job (not recommended), you must create an encrypted MycotoolsDB passkey using `updateDB.py` and pass the password to stdin to these scripts.

<br />

These scripts can input assembly accessions (NCBI) or genome codes (JGI). The column must have the appropriate header ('assembly_acc' or 'assembly_acc') (substitute `-i aspergillus.mtdb` with this file):

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

You can download NCBI SRA's by acquiring NCBI's SRA tools, making sure `fastq-dump` is included in your PATH, and then running:

```bash
ncbiDwnld.py --sra -i <REFERENCE>
```

You can create a file with SRA ID's or BioProject, etc. Basically any query that is unique and sufficient to acquire the SRRs of interest. For paired-end reads, append `-pe` to the command.

<br /><br />


## Grab accessions
### acc2fa.py / acc2gff.py / acc2gbk.py
By default, if you are querying using a MycotoolsDB accession then it can search the database without a standalone file.
Let's say you want to query *Panaeolus cyanescens'* PsiD and the NCBI accession is "PPQ80975.1". Find Panaelous cyanescens' ome code in the database:
```bash
grep Panaeolus $(mycodb) | grep cyanescens | cut -f 1
```

The first column in the output is `pancya1`, which is the ome code for this organism. Now, prepend the code to the accession and grab it from the db (you can use `>` after the commands to pipe output to a file):
```bash
acc2gff.py -a pancya1_PPQ80975.1
acc2fa.py -a pancya1_PPQ80975.1
acc2gbk.py -a pancya1_PPQ80975.1
```

If you have a list of accessions, create an input file with the accessions separated by new lines then run:
```bash
acc2gff.py -i <INPUTFILE>
acc2fa.py -i <INPUTFILE>
acc2gbk.py -i <INPUTFILE>
```

<br /><br />
## Fasta coordinates
### coords2fa.py

```bash
coords2fa.py -h

Input nucleotide fasta/tsv input, extract coordinates
coords2fa.py <FA> <SEQID> <START_COORD> <END_COORD> <STRAND_SENSE>
Extract full sequence from sense strand: coords2fa.py test.fna scaffold_20 0 -1
Extract coordinates from antisense strand: coords2fa.py test.fna scaffold_20 69 420 -


Bulk extraction tab delimitted row format:
fasta   sequence_id     start_coordinate        end_coordinate  strand_sense

coords2fa.py coords.tsv
```

<br /><br />

## Gene coordinates to sequences
### gff2seq.py
`gff2seq.py` will extract the nucleotide or amino acid sequences associated with a gene coordinates `gff` file. Users can optionally input a flanking nucleotide plus/minus to extract from 
the flanks of the provided `gff`. `gff2seq.py` will extract flanks independently for each sequence ID (column 1) within the `gff`. Coding or noncoding regions can be specified. If `--intergenic` is called then only the first and last gene for each fasta
sequence are considered.

e.g. extract nucleotide sequences and 1 kilobase flanks and noncoding regions within the following genes:

```bash
gff2seq.py -g <.GFF3> -a <.FNA> -nc -pm 1000 -n
```

```
optional arguments:
  -h, --help            show this help message and exit
  -g GFF, --gff GFF
  -n, --nucleotide
  -p, --protein
  -a ASSEMBLY, --assembly ASSEMBLY
  -i --intergenic	-n only
  -nc, --noncoding      -n only
  -pm PLUSMINUS, --plusminus PLUSMINUS
                        -n only
  -af, --all_flanks     -n and -nc only
```

<br /><br />


## Grab loci
### acc2loci.py
Inputs an accession (`-a`) or new line separated list of accessions and
optional genes +/- (`-p`, default 10). Outputs gene accessions or a gff and
protein fasta of the loci.

output gff and protein fasta of an accession's cluster (outputs to `<ACCESSION>_clus*`):
```bash
acc2loci.py -o -a <MYCOTOOLSDB_ACCESSION>
```

list proximal +/- 5 genes to standard out:
```bash
acc2loci.py -a fibsp.1_906341 -p 5

fibsp.1_880711
fibsp.1_846234
fibsp.1_809145
fibsp.1_923701
fibsp.1_771516
fibsp.1_906341
fibsp.1_719531
fibsp.1_846242
fibsp.1_138
fibsp.1_942299
fibsp.1_906343
```

<br /><br />

## Curate annotation
### curAnnotation.py
`curAnnotation.py` is tailored toward curating OrthoFiller or Funannotate output (more curation available upon request). This script will convert OrthoFiller `.gtf` to `.gff3`, rename headers sequentially, and add an `Alias=<PREFIX>` field for MycotoolsDB compatible accession for each entry.

```bash
curAnnotation.py -g <ORTHOFILLER>/results/results.gtf -f <ASSEMBLY> -p <CODENAME>
```

### curGFF3.py / curProteome.py / gff2gff3.py
There are several scripts in the `utils` used to curate gene coordinate files and proteomes for the MycotoolsDB. `curGFF3.py` is tested with both JGI and NCBI `gff3` files, `gff2gff3.py` curates JGI `gff2` files to MycotoolsDB compatible `gff3`, and `curProteome.py` curates NCBI or JGI proteomes.

<br /><br />

## Visualizing Loci
### gff2svg.py
This script will input a .gff3, or new line-delimited list of .gff3 paths, 
and for each contig output an SVG of the locus annotated by function retrieved 
from the 'product=' field. 

<br />

e.g. make an SVG from a GFF3: `gff2svg.py -g <MY.gff3>`

make SVGs for all GFF3s in a new line delimited list with width set to 20:

```bash
gff2svg.py -i <LISTOFGFF3.nsv> -o <OUTPUT_DIR> -w 20
```


<br /><br /><br />

# EVOLUTIONARY ANALYSIS TOOLS
## BLAST MycotoolsDB
### db2search.py
`db2search.py` will `blastn`, `blastp`, `tblastn`, or `blastx` the MycotoolsDB using a query fasta and compile a results fasta for each accession in the query according to any inputted threshold requirements. It is recommended to keep `--cpu` below the number the number of query organisms.

```bash
db2search.py --help
usage: db2search.py [-h] -b BLAST [-q QUERY] [-e EVALUE] [-s BITSCORE] [-i IDENTITY] [-m MAXHITS] [-d DATABASE]
                   [-o OUTPUT] [-p PREVIOUS] [--cpu CPU]

Blasts a query against a db and compiles output fastas for each query.

optional arguments:
  -h, --help            show this help message and exit
  -b BLAST, --blast BLAST
                        Blast type { blastn, blastp, tblastn, blastx }
  -q QUERY, --query QUERY
                        Query fasta
  -e EVALUE, --evalue EVALUE
                        Negative e-value order max, e.g. 2 for 10^-2
  -s BITSCORE, --bitscore BITSCORE
                        Bit score minimum
  -i IDENTITY, --identity IDENTITY
                        Identity minimum, e.g. 0.6
  -m MAXHITS, --maxhits MAXHITS
                        Max hits for each organism
  -d DATABASE, --database DATABASE
                        mycotoolsdb, DEFAULT: masterdb
  -o OUTPUT, --output OUTPUT
  -p PREVIOUS, --previous PREVIOUS
                        Previous run directory
  --cpu CPU             DEFAULT: all
```

<br /><br />


## hmmsearch MycotoolsDB
### db2hmm.py
`db2hmm.py` will compile hmmsearch results, optionally ouput fasta/hmmalign to original models/trim alignments from a profile hidden markov model. This script supports multiprocessing and uses all detected cores by default.

e.g.: `db2hmm.py -d profile.hmm` will run the master database referencing the inputted hmmdb. 

`db2hmm.py -d profile.hmm -b 1 -t 75` takes the best hit from hits with 25 - 100% query coverage
```bash
db2hmm.py --help
Runs `hmmsearch` v each proteome in a `.mtdb`. For each query, extracts results
and optionally outputs a compiled fasta, hmmalignment and/or trimmed
alignment.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input database `.mtdb`
  -d HMMDB, --hmmdb HMMDB
                        Hmm database `.hmm`
  -o OUTPUT, --output OUTPUT
                        User-specified output directory
  -p PREVIOUS, --previous PREVIOUS
                        Previous db2hmm dir. Be wary of incomplete
                        outputs from interrupted runs
  -c CPU, --cpu CPU     Processors to use. Default = All
  -f, --fasta           Compile fastas for each query
  -l, --align           Align fastas to original hmm and trim via `trimal`.
                        Calls `-f`
  -a, --accession       Extract accessions instead of queries (Pfam, etc).
                        Requires `-f`
  -b BEST, --best BEST  # top hits for each organism Requires `-f`
  -t THRESHOLD, --threshold THRESHOLD
                        Query percent threshold (+/-). Requires `-f`
  -e EVALUE, --evalue EVALUE
                        E value threshold, e.g. 10^(-x) where x is the input.
                        Requires `-f`
  --trimal TRIMAL       User-specified trimAl commands in "", e.g.:
                        "-strictplus -fasta"
```

<br /><br />


## Tree Building
### fa2tree.py
`fa2tree.py` will input a fasta file or directory of fasta files, trim, and
generate trees either via job submission (`slurm` or `torque`) or immediate execution. 
Note you will need `mafft`, `clipkit`, and `iqtree` installed. 
If these are not installed, install them into your mycotools conda environment via 
`conda install mafft iqtree` and `pip install clipkit`.

Information on a complete [phylogenetic pipeline](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#mycotools-pipelines) is elaborated below.

<br />

e.g. Swiftly construct a tree from a fasta
```bash
fa2tree.py -f <FASTA>.fa --fast
```

Prepare a robust tree for slurm job submission
```bash
fa2tree.py -f <FASTA> -s -A PAS1046
```

Begin the tree pipeline by navigating into the folder and running `bash
mafft.sh` (execute immediately) or `sbatch mafft.sh` (job submission). If you
run out of memory, then add increased memory/cores to the job submission script(s).

Most often, the final file you want is `*.contree`.
View your trees using [FigTree](https://github.com/rambaut/figtree/releases).

<br /><br />


## Hierarchical agglomerative clustering
### fa2clus.py
Hierarchical agglomerative clustering is a useful systematic approach to
clustering groups of sequences for phylogenetic analysis via percent identity.
Some gene families (e.g. P450s) yield 10,000s of results for BLAST searches against the
MycotoolsDB. `fa2clus.py` allows the user to truncate these sets without constructing a
phylogeny using a minimum percent identity cutoff. `fa2clus.py` optionally implements
an automated iterative approach to obtaining a cluster of minimum - maximum size
with the gene of interest.

`fa2clus.py` will either take a `fasta` and generate a distance matrix using 
`usearch calc_distmx` by default or the % identity of `needle` alignments.
Then, cluster sequences via hierarchical agglomerative clustering and output a
`.clus` file of cluster assignments and `.newick` dendrogram. 

Currently, using `needle` is slow, so it is recommended to acquire a
free license for `usearch` and use that instead. The limitations of the free
license are sufficient for distance matrix calculation even on large datasets. 

e.g. Calculate a distance matrix and cluster from a fasta with a minimum
identity 0.3 and maximum distance 0.7 (1 - identity) to consider a connection:
```bash
fa2clus.py -f <FASTA>.fa -m 0.2 -x 0.7
```

Iteratively cluster until a cluster size of 50-200 genes is achieved:
```bash
fa2clus.py -f <FASTA> -m 0.2 -x 0.7 --iterative <GENE> --minseq 50 --maxseq 200
```


<br /><br /><br />


# Mycotools Pipelines
## Phylogenetic analysis

A key component of robust phylogenetic reconstruction is
adequately sampling. MycotoolsDB enables adequate sampling by 
providing a near-comprehensive database of available fungal genomic data. 
Despite the benefits of increased sampling, there are two prominent problems: 1) large alignments
lose resolution because they must consider many sequences and 2) computational 
complexity increases with sample size. 
For a small set of genes, such as ITS, one can usually assume that the gene family is 
strictly vertically inherited, so the dataset can be cut down to closely 
related organisms. For most other genes, it is not valid to assume the gene family
is vertically conserved because horizontal transfer is a prominent modality of
gene transmission in prokaryotes and fungi.

It is thus both important to adequately sample and systematically truncate 
the dataset into a manageable set of gene homologs. This is accomplished by
iteratively constructing phylogenies, identifying gene family homologs,
truncating the data to these homologs, and repeating until a manageable tree is
obtained. On its own, this analysis requires elaborate integration of multiple 
independent softwares, but Mycotools takes care of the bulk of this work.

You will need several programs for this analysis, so activate your Mycotools conda environment
and run `conda install iqtree mafft blast -c bioconda` and `pip install clipkit`

<br />

### Example 1:

A phylogenetic analysis often starts with a single gene of interest. The first step
is to obtain the gene family of closely related genes by BLASTing a gene protein sequence 
across the mycotoolsDB. 

1. extract a database of published sequences, or use other arguments to extract
other organisms of interest

```bash
extractDB.py > pub.mtdb
```

<br />

2. obtain gene homologs using `db2search.py` with an e-value threshold of 10<sup>-2</sup>:

```bash
db2search.py -q <QUERYGENE.fasta> -b blastp -e 2 -d pub.mtdb
```

<br />

3. `db2search.py` created a directory with fasta results. obtain the number of
genes recovered from the analysis:

```bash
grep '>' <RESULTS.fasta> | wc -l
```

<br />

4a. If there are fewer than 10,000 sequences you can proceed directly to robust tree
building. Otherwise, proceed to the `b` steps. If the sample size is too large
to finish the analysis, specify `--fast` below. 

```bash
fa2tree.py -f <BLASTRESULTS.fasta> -s -A
<PROJECT>
```

<br />

5a. `fa2tree.py` created a directory, now submit the job

```bash
sbatch <TREEDIRECTORY>/mafft.sh
```

<br />

6a. Download the tree (usually `*.contree` for iqtree). Open
the tree in [figtree](https://github.com/rambaut/figtree/releases). If the
dataset is too big, identify a highly supported node (no less than 95, ideally
100), containing your sequences of interest. Change to tip mode at the top of
the program, click on the node, and copy the sequences highlighted into a plain
text file. Then, restart the phylogenetic reconstruction by obtaining the
sequences using `acc2fa.py` and restarting 4a with the output fasta.

```bash
acc2fa.py -i <FILEWITHACCESSIONS>
```

<br />

4b. If there are greater than 10,000 sequences, trees - even fasttrees can be
intractable at this quantity. Therefore, try hierarchical agglomerative
clustering. Adjust the % identity for both arguments as needed. You may need to
submit this as a job

```bash
fa2clus.py -f <FASTA>.fa -m 0.3 -x 0.7
```

<br />

5b. Open the `.clus` file to obtain the cluster with your sequence of interest,
copy those sequences into a blank file, then use `acc2fa.py` as described in
step 6a. If the clusters are too small/big, open the `.newick` output in figtree to 
select the clade of interest as described in step 6a.

<br />

6b. Repeat the analysis at step 4a until a final phylogeny is obtained.


<br /><br />


## Cluster Reconstrunction and Phylogenetic Analysis (CRAP)
### crap.py
CRAP (originally created by Jason Slot) is a simple, elegant pipeline for
studying the evolution of a gene cluster on a gene-by-gene basis. CRAP 
will 1) input a cluster query and use a search algorithm (BLAST/mmseqs/Diamond)
or orthogroup-based approach to find homologs in the MycotoolsDB; 2) use 
hierarchical clustering to truncate the sequence set and detect outgroups;
3) construct phylogenies of each query sequence; and 4) map locus synteny 
diagrams  onto the leaves of the phylogenies.

`crap.py` can operate on a query of MycotoolsDB accessions or a standalone
multifasta input of external accessions. Following homolog acquisition,
`crap.py` will submit each set of hits for tree building or hierarchical 
agglomerative clustering if the number of sequences exceeds the inputted maximum 
(`-m`).

By default, `crap.py` will construct trees using fasttree. Fasttree is usually
sufficient to get an idea of how the cluster is evolving because CRAP builds
phylogenies for all query genes within a cluster, so congruent topology across
different genes is a good window into the evolution. Alternatively, CRAP can
also construct 1000 bootstrap iterations
of IQTree with the ModelFinder module by specifying `-i`/`--iqtree`.

To search an extracted sub-MycotoolsDB using `blastp` and create phylogenies with `fasttree`:
```bash
extractDB.py --rank phylum --lineage Basidiomycota > basi.mtdb
crap.py -q <QUERYGENES> -d basi.mtdb -s blastp --bitscore 40 --cpu 12
```

<br />

```bash
usage: crap.py [-h] -q QUERY [-d DATABASE] [-s SEARCH] [-b BITSCORE] [-og ORTHOGROUPS] [-p PLUSMINUS]
               [--maxseq MAXSEQ] [-i] [--conversion CONVERSION] [--ingroup] [--no_label] [-g GFF]
               [-o OUTPUT] [-c CPU] [-v]

Mycotools integrated Cluster Reconstruction and Phylogeny (CRAP) pipeline

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY, --query QUERY
                        Fasta or white-space delimited file/string of cluster genes
  -d DATABASE, --database DATABASE
  -s SEARCH, --search SEARCH
                        Search binary {mmseqs, blastp} for search-based CRAP
  -b BITSCORE, --bitscore BITSCORE
                        Bitscore minimum for search algorithm. DEFAULT: 30
  -og ORTHOGROUPS, --orthogroups ORTHOGROUPS
                        MycotoolsDB Orthogroup tag for OG-based CRAP. DEFAULT: P for phylum
  -p PLUSMINUS, --plusminus PLUSMINUS
                        Genes up-/downstream to analyze from loci. DEFAULT: 10
  --maxseq MAXSEQ       Max sequences for trees/min for fa2clus.py. Outgroup detection may exceed this
                        number. DEFAULT: 250
  -i, --iqtree          IQTree2 1000 bootstrap iterations. DEFAULT: fasttree
  --conversion CONVERSION
                        Tab delimited conversion file for annotations: Query Conversion
  --ingroup             Do not detect outgroups, do not root trees
  --no_label            Do not label synteny diagrams
  -g GFF, --gff GFF     GFF for non-mycotools locus diagram. Requires -s and a fasta for -i
  -o OUTPUT, --output OUTPUT
                        Output base dir. Will rerun if previous director exists.
  -c CPU, --cpu CPU
  -v, --verbose
```

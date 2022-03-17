# INDEX

<br />


- **MYCOTOOLSDB TOOLS**
	- [Initializing MycotoolsDB](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#initialization)
	- [Updating MycotoolsDB](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#updating]
	- [Interfacing with the database](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#mycotools-db)
	- [Creating modular databases](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#creating-modular-databases)
	- [Acquiring database files / file paths](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#acquiring-database-files)
	- [Substitute organism name for MycotoolsDB organism code](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#ome2name.py)

<br />


- **SEQUENCE DATA**
	- [Downloading from NCBI / JGI](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#downloading-files)
	- [Sequence data statistics](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#sequence-data-statistics)
	- [Grabbing accessions](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#grab-accessions)
	- [Grabbing loci](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#grab-clusters)
	- [Visualizing loci](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#visualizing-loci)
	- [Curating annotation](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#curate-annotation)

<br />


- **EVOLUTIONARY ANALYSES**
	- [MycotoolsDB BLAST](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#blast-mycotoolsdb)
	- [MycotoolsDB hidden markov model search](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#hmmsearch-mycoDB)
	- [Fasta to tree](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#tree-building)
	- [Hiearchical agglomerative clustering](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#hiearchical-agglomerative-clustering)

<br />

- **MYCOTOOLS PIPELINES**
    - [Phylogenetic analysis](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#mycotools-pipelines)


---

<br /><br /><br />


# MYCOTOOLS DB
The following are information regarding scripts that manipulate Mycotools .db files. To learn more about Mycotools .db, you may refer to the master [MycotoolsDB repository](https://gitlab.com/xonq/mycotoolsdb/-/README.md) (currently unavailable).

<br /><br />

## Initialization
### updateDB.py
To initialize a curated database of all NCBI and MycoCosm (JGI) genomes, run the following

```bash
(mycotools) -$ updateDB.py -i <INIT DIRECTORY>
```

<br /><br />

## Updating
To update MycotoolsDB:

```bash
(mycotools) -$ updateDB.py -u
```

<br /><br />

## Interfacing with the master database
### mycodb
`mycodb` is a utility that integrates with the master database or just prints the path of the master database, which can then be used with other shell commands. MycotoolsDBs are labelled `YYYYmmdd.db`.
```bash
(mycotools) -$ mycodb
/home/xonq/mycodb/mycodb/20210125.db
```

If you want to use the path of the master database you can use basic bash functionality to work with the output, 
e.g. to open in a text editor or to grep the file:
```bash
(mycotools) -$ vim $(mycodb)
(mycotools) -$ grep 'Psilocybe' $(mycodb)
```

<br /><br />


## Creating modular databases
### extractDB.py
If you are only interested in a subset of lineages in the master mycotoolsDB, then extract the portion you want, run `extractDB.py -h` to see all options:

e.g. grab a database of a taxonomic order: 
```bash
(mycotools) -$ extractDB.py -l Atheliales -r order > atheliales.db
```

grab all NCBI Aspergilli accessions: 
```bash
(mycotools) -$ extractDB.py -s ncbi -l aspergillus -r genus > aspergillus.db_ncbi
``` 

grab the inverse of your arguments: 
```bash
(mycotools) -$ extractDB.py -s ncbi -l aspergillus -r genus -i > notAspergullis.db_notNcbi
```

grab a list of orders from a file:
```bash
(mycotools) -$ extractDB.py -ll <TAX_FILE> -r order > taxa.db
```

grab a list of `ome`s in a new line delimited file: 
```bash
(mycotools) -$ extractDB.py ---ome <OME_FILE>
```

```bash
(mycotools) -$ extractDB.py --help
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
### dbFiles.py
Inputs a MycotoolsDB `.db` file (by default uses the master database), then creates symlinks of the selected file types, hard copies the files, or prints their PATHs. A symlink is simply creating a placeholder file that links to the database file... this way it does not take up additional storage space like a hard copy does. However, editing symlinks will edit the original file, so *only hard copy `--hard` if you need to edit the files*.

Let's say you want protein data from organisms in one family. First, you should extract a database of organisms you want:
```bash
(mycotools) -$ mkdir pullFiles && cd pullFiles
(mycotools) -$ extractDB.py -r family -l Atheliaceae > atheliaceae.db
```

Then, run `dbFiles.py` to copy the protein fastas into the current directory (call `-h` to see all options):
```bash
(mycotools) -$ dbFiles.py -d atheliaceae.db -p 
```

Alternatively, if you just need the paths (links) to these files, simply run:
```bash
(mycotools) -$ dbFiles.py -d atheliaceae.db -p --print
```

<br /><br />


## Other MycotoolsDB scripts
### ome2name.py
Substitutes MycotoolsDB organism code names (e.g. `fusgra1`) for taxonomic information (e.g. Fusarium_graminearum_var._XYZ).

e.g. to substitute ome for genus species and strain: `ome2name.py <INPUT> oa`
```bash
(mycotools) -$ ome2name.py -h
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
(mycotools) -$ assemblyStats.py <ASSEMBLY.fa>
(mycotools) -$ annotationStats.py <ANNOTATION.gff3>
```

To obtain a table of organisms' annotation statistics, [create a mycotoolsDB](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#creating-modular-databases) file with the organisms of interest and run:
```bash
(mycotools) -$ assemblyStats.py <MYCOTOOLSDB.db>
(mycotools) -$ annotationStats.py <MYCOTOOLSDB.db>
```

If you want to route the output to a file, simply redirect output by appending ` > <OUTPUTFILE>` to the command, or add an output file as the second argument


<br /><br />


## Downloading files
### jgiDwnld.py / ncbiDwnld.py
These scripts input a MycotoolsDB or can be manually made as shown at the bottom of this section. 

Say you want to grab a few organisms' transcript information from your genus, *Aspergillus*. First, extract entries in the database that are within *Aspergillus*:
```bash
(mycotools) -$ mkdir dwnldFiles && cd dwnldFiles
(mycotools) -$ extractDB.py -c genus -t aspergillus > aspergillus.db_ncbi
```

If there are organisms you don't want in the extracted `.db`s, just delete their line(s) in the file. Next call `jgiDwnld.py -h` or `ncbiDwnld.py -h` to find the flags necessary to download the files you want. To download transcript data (and EST data for JGI) in your current directory:
```bash
(mycotools) -$ jgiDwnld.py -i aspergillus.db_jgi -t -e
(mycotools) -$ ncbiDwnld.py -i aspergillus.db_ncbi -t
```

These scripts populate with compressed files. To unzip all the files, run `gunzip <FILETYPE>/*.gz`. You will also see log files for the download process. 
To submit as a job (not recommended), you must create an encrypted MycotoolsDB passkey using `updateDB.py` and pass the password to stdin to these scripts.

<br />

These scripts can input biosamples (NCBI) or genome codes (JGI). The column must have the appropriate header ('genome_code' or 'biosample') (substitute `-i aspergillus.db` with this file):

`jgiGenomeCodes.txt`
```
genome_code
Abobi1
Absrep1
Acain1
```

`ncbiBioSamples.txt`
```
biosample
SAMN15352002
SAMN06205226
SAMN02744098
```

<br />

You can download NCBI SRA's by acquiring NCBI's SRA tools, making sure `fastq-dump` is included in your PATH, and then running:

```bash
(mycotools) -$ ncbiDwnld.py --sra -i <REFERENCE>
```

You can create a file with SRA ID's or BioProject, etc. Basically any query that is unique and sufficient to acquire the SRRs of interest. For paired-end reads, append `-pe` to the command.

<br /><br />


## Grab accessions
### acc2fa.py / acc2gff.py
By default, if you are querying using a MycotoolsDB accession then it can search the database without a standalone file.
Let's say you want to query *Panaeolus cyanescens'* PsiD and the NCBI accession is "PPQ80975.1". Find Panaelous cyanescens' ome code in the database:
```bash
(mycotools) -$ grep Panaeolus $(mycodb) | grep cyanescens | cut -f 1
```

The first column in the output is `pancya1`, which is the ome code for this organism. Now, prepend the code to the accession and grab it from the db (you can use `>` after the commands to pipe output to a file):
```bash
(mycotools) -$ acc2gff.py -a pancya1_PPQ80975.1
(mycotools) -$ acc2fa.py -a pancya1_PPQ80975.1
```

If you have a list of accessions, create an input file with the accessions separated by new lines then run:
```bash
(mycotools) -$ acc2gff.py -i <INPUTFILE>
(mycotools) -$ acc2fa.py -i <INPUTFILE>
```

<br /><br />


## Grab loci
### grabLoci.py
Inputs an accession (`-a`) or new line separated list of accessions and optional genes +/- (`-p`, default 10). Outputs gene accessions or a gff and protein fasta of the clusters.

output gff and protein fasta of an accession's cluster (outputs to `<ACCESSION>_clus*`):
```bash
(mycotools) -$ grabLoci.py -o -a <MYCOTOOLSDB_ACCESSION>
```

list proximal +/- 5 genes to standard out:
```bash
(mycotools) -$ grabLoci.py -a fibsp.1_906341 -p 5

fibsp.1_906341 locus +/- 5
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
(mycotools) -$ curAnnotation.py -g <ORTHOFILLER>/results/results.gtf -f <ASSEMBLY> -p <CODENAME>
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
(mycotools) -$ gff2svg.py -i <LISTOFGFF3.nsv> -o <OUTPUT_DIR> -w 20
```


<br /><br /><br />

# EVOLUTIONARY ANALYSIS TOOLS
## BLAST MycotoolsDB
### db2blast.py
`db2blast.py` will `blastn`, `blastp`, `tblastn`, or `blastx` the MycotoolsDB using a query fasta and compile a results fasta for each accession in the query according to any inputted threshold requirements. It is recommended to keep `--cpu` below the number the number of query organisms.

```bash
(mycotools) -$ db2blast.py --help
usage: db2blast.py [-h] -b BLAST [-q QUERY] [-e EVALUE] [-s BITSCORE] [-i IDENTITY] [-m MAXHITS] [-d DATABASE]
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
### db2hmmsearch.py
`db2hmmsearch.py` will compile hmmsearch results, optionally ouput fasta/hmmalign to original models/trim alignments from a profile hidden markov model. This script supports multiprocessing and uses all detected cores by default.

e.g.: `db2hmmsearch.py -d profile.hmm` will run the master database referencing the inputted hmmdb. 

`db2hmmsearch.py -d profile.hmm -b 1 -t 75` takes the best hit from hits with 25 - 100% query coverage
```bash
(mycotools) -$ db2hmmsearch.py --help
Runs `hmmsearch` v each proteome in a `.db`. For each query, extracts results
and optionally outputs a compiled fasta, hmmalignment and/or trimmed
alignment.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input database `.db`
  -d HMMDB, --hmmdb HMMDB
                        Hmm database `.hmm`
  -o OUTPUT, --output OUTPUT
                        User-specified output directory
  -p PREVIOUS, --previous PREVIOUS
                        Previous db2hmmsearch dir. Be wary of incomplete
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
Note you will need to have `mafft`, `trimal`, and the tree software you are
using (`fasttree`/`iqtree`) installed. If these are not installed, install them
into your mycotools conda environment via `conda install mafft trimal iqtree`

e.g. Construct a tree from a fasta
```bash
(mycotools) -$ fa2tree.py -i <FASTA>.fa -t fasttree
```

Prepare for slurm job submission
```bash
(mycotools) -$ fa2tree.py -i <FASTA> -t iqtree -s -A PAS1046
```
https://github.com/rambaut/figtree/releases
Begin the tree pipeline by navigating into the folder and running `bash
mafft.sh` (execute immediately) or `sbatch mafft.sh` (job submission).

<br /><br />


## Hierarchical agglomerative clustering
### aggClus.py
Hierarchical agglomerative clustering is a useful systematic approach to
extracting groups of sequences for phylogenetic analysis. Ubiquitous genes,
like P450s, will often yield 10,000s of results for BLAST searches against the
MycotoolsDB. Constructing and visualizing a tree of this magnitude just is not
practical in many cases; it is therefore necessary to decrease the size of the 
dataset to a workable size while relying on biological information (global
pairwise alignments).

`aggClus.py` will either take a `fasta` and generate a distance matrix using 
`usearch calc_distmx` by default or the % identity of `needle` alignments.
Then, cluster sequences via hierarchical agglomerative clustering and output a
`.clus` file of cluster assignments and `.newick` dendrogram. 

Currently, using `needle` takes quite long, so it is recommended to acquire a
free license for `usearch` and use that instead. The limitations of the free
license are sufficient for distance matrix calculation even on large datasets. 

e.g. Calculate a distance matrix and cluster from a fasta with a minimum
identity 0.3 and maximum distance 0.7 (1 - identity) to consider a connection:
```bash
(mycotools) -$ aggClus.py -f <FASTA>.fa -m 0.3 -x 0.7
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
is vertically conserved because horizontal transfer is a prominent modality of gene 
evolution in Fungi. 

It is thus both important to adequately sample and systematically truncate 
the dataset into a manageable set of gene homologs. This is accomplished by
iteratively constructing phylogenies, identifying gene family homologs,
truncating the data to these homologs, and repeating until a manageable tree is
obtained. On its own, this analysis requires elaborate integration of multiple 
independent softwares, but Mycotools takes care of the bulk of this work.

You will need several programs for this analysis, so activate your Mycotools conda environment
and run `conda install iqtree mafft trimal blast -c bioconda` 

<br />

### Example 1:

A phylogenetic analysis often starts with a single gene of interest. The first step
is to obtain the gene family of closely related genes by BLASTing a gene protein sequence 
across the mycotoolsDB. 

1. extract a database of published sequences, or use other arguments to extract
other organisms of interest

```bash
extractDB.py > pubFungi.db
```

<br />

2. obtain gene homologs using `db2blast.py` with an e-value threshold of 10<sup>-2</sup>:

```bash
(mycotools) -$ db2blast.py -q <QUERYGENE.fasta> -b blastp -e 2 -d pubFungi.db
```

<br />

3. `db2blast.py` created a directory with fasta results. obtain the number of
genes recovered from the analysis:

```bash
(mycotools -$ grep '>' <RESULTS.fasta> | wc -l
```

<br />

4a. If there are fewer than 10,000 sequences you can proceed directly to tree
building. Otherwise, proceed to the `b` steps. If the sample size is too large
to finish the analysis, specify `fasttree` below; however, if this is the final 
tree, or if the dataset is small enough, use `iqtree` because it is more robust:

```bash
(mycotools) -$ fa2tree.py -f <BLASTRESULTS.fasta> -t <TREESOFTWARE> -s -A
<PROJECT>
```

<br />

5a. `fa2tree.py` created a directory, now submit the job

```bash
(mycotools) -$ sbatch <TREEDIRECTORY>/mafft.sh
```

<br />

6a. Download the tree (`*.contree` for iqtree, `*.fasttre` for fasttree). Open
the tree in [figtree](https://github.com/rambaut/figtree/releases). If the
dataset is too big, identify a highly supported node (no less than 95, ideally
100), containing your sequences of interest. Change to tip mode at the top of
the program, click on the node, and copy the sequences highlighted into a plain
text file. Then, restart the phylogenetic reconstruction by obtaining the
sequences using `acc2fa.py` and restarting 4a with the output fasta.

```bash
(mycotools) -$ acc2fa.py -i <FILEWITHACCESSIONS>
```

<br />

4b. If there are greater than 10,000 sequences, trees - even fasttrees can be
intractable at this quantity. Therefore, try hierarchical agglomerative
clustering. Adjust the % identity for both arguments as needed. You may need to
submit this as a job

```bash
(mycotools) -$ aggClus.py -f <FASTA>.fa -m 0.3 -x 0.7
```

<br />

5b. Open the `.clus` file to obtain the cluster with your sequence of interest,
copy those sequences into a blank file, then use `acc2fa.py` as described in
step 6a. If the clusters are too small/big, open the `.newick` output in figtree to 
select the clade of interest as described in step 6a.

<br />

6b. Repeat the analysis at step 4a until a final phylogeny is obtained.

# INDEX

<br />

- **MYCODB**
	- [Creating modular databases](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#creating-modular-databases)
	- [Acquiring database files / file paths](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#acquiring-database-files)

<br />

- **SEQUENCE DATA**
	- [Downloading from NCBI / JGI](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#downloading-files)
	- [Sequence data statistics](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#sequence-data-statistics)
	- [Grabbing accessions](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#grab-accessions)

<br />

- **ANALYSES**
	- [mycoDB blast](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#blast-mycodb)
	- [mycoDB hidden markov model search](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#hmmsearch-mycoDB)
	- [fasta to tree](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#phylogenetic-analysis)
	- [Hiearchical agglomerative clustering](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md#hierarchical-agglomerative-clustering)

---
<br /><br />

## Creating modular databases
### abstractDB.py
The core principle behind mycotools is modularity - creating tools that modularly interface with databases ranging from 1 organism to all ~ 1,500 published fungi. Most scripts use the master database by default, but if you are only interested in a portion of the database, then abstract the portion you want.

e.g. grab a database of a taxonomic order you're interested in: `abstractDB.py -t Atheliales -c order > atheliales.db`

grab all NCBI Aspergilli accessions: `abstractDB.py -s ncbi -t aspergillus -c genus > aspergillus.db_ncbi` 

grab a list of orders in a new line delimited file: `abstractDB.py -tl <TAX_FILE> -c order > taxa.db`
```
Imports database and file of omes or taxonomy. Abstracts database from
parameters. E.g. `abstractDB.py -i $DATABASE -t Atheliaceae -c family`

optional arguments:
  -h, --help            show this help message and exit
  -t TAXONOMY, --taxonomy TAXONOMY
                        Taxonomic group to abstract. Requires `-c`
  -c CLASSIFICATION, --classification CLASSIFICATION
                        Taxonomic classification to abstract from. Requires
                        `-t` or `-tl`
  -s SOURCE, --source SOURCE
                        Genome source to abstract
  -n, --nonpublished    Abstract nonpublished rows.
  --ome OME             New line separated list of omes to include.
  -tl TAXONOMY_LIST, --taxonomy_list TAXONOMY_LIST
                        New line separated list of taxonomic groups - must be
                        same classification. Requires `-c`
  --inverse             Inverse arguments
  --headers             Do not output db headers
  -i INPUT, --input INPUT
                        Kontools database. DEFAULT: master database
  -o OUTPUT, --output OUTPUT
                        Output path instead of print to stdout. Includes
                        column headers (stdout does not)
```

<br />

## Acquiring database files
### dbFiles.py
Inputs a mycotools `.db` file (by default uses the master database), then pulls the selected file types or prints their PATHs.

Let's say you want protein data from organisms in one family. First, you should abstract a database of organisms you want:
```
mkdir pullFiles && cd pullFiles
abstractDB.py -c family -t Atheliaceae > atheliaceae.db
```

Then, run `dbFiles.py` to copy the protein fastas into the current directory (call `-h` to see all options):
```
dbFiles.py -i atheliaceae.db -p 
```

Alternatively, if you just need the paths (links) to these files, simply run:
```
dbFiles.py -i atheliaceae.db -p -l
```

<br /><br />

## Sequence data statistics
### assemblyStats.py / annotationStats.py
```
assemblyStats.py <ASSEMBLY.fa>
annotationStats.py <ANNOTATION.gff3>
```

<br />

## Downloading files
### jgiDwnld.py / ncbiDwnld.py
These scripts` input can be obtained from the mycotools database as described below, or can be manually made as shown at the bottom of this section. Say you want to grab a few organisms' transcript information from your genus, *Aspergillus*. First, abstract entries in the database that are within *Aspergillus* for both JGI & NCBI:
```
mkdir dwnldFiles && cd dwnldFiles
abstractDB.py -s jgi -c genus -t aspergillus > aspergillus.db_jgi
abstractDB.py -s ncbi -c genus -t aspergillus > aspergillus.db_ncbi
```

If there are organisms you don't want in the abstracted `.db`s, just delete their line(s) in the file. Next call `jgiDwnld.py -h` or `ncbiDwnld.py -h` to find the flags necessary to download the files you want. To download transcript data (and EST data for JGI) in your current directory run:
```
jgiDwnld.py -i aspergillus.db_jgi -t -e
ncbiDwnld.py -i aspergillus.db_ncbi -t
```

You will now see folders named after the file types you downloaded and the compressed files stored within. To unzip all the files, run `gunzip <FILETYPE>/*.gz`. You will also see log files for the download process. Don't submit these as a job, you are required to enter a password for `jgiDwnld.py`, and if the command stops you can simply rerun in the same directory and it should take-off where you left it/where it failed.

<br />

Manually created input files - the important thing is the column has the exact appropriate header ('genome_code' or 'basename') (substitute `-i aspergillus.db` with this file):

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

## Grab accessions
### acc2fa.py / acc2gff.py
By default, if you are querying using a mycotools accession then it can search the database without a standalone file.
Let's say you want to query *Panaeolus cyanescens'* PsiD and the NCBI accession is "PPQ80975.1". Find Panaelous cyanescens' ome code in the database:
```
grep Panaeolus $(mycodb) | grep cyanescens | cut -f 1
```

The first column in the output is `pancya1`, which is the ome code for this organism. Now, prepend the code to the accession and grab it from the db (you can use `>` after the commands to pipe output to a file):
```
acc2gff.py -a pancya1_PPQ80975.1
acc2fa.py -a pancya1_PPQ80975.1
```

If you have a list of accessions, create an input file with the accessions separated by new lines then run:
```
acc2gff.py -i <INPUTFILE>
acc2fa.py -i <INPUTFILE>
```


<br /><br />

## BLAST mycoDB
### db2blast.py
`db2blast.py` will `blastn`, `blastp`, `tblastn`, or `blastx` the mycoDB using a query fasta and compile a results fasta for each accession in the query according to any inputted threshold requirements. It is recommended to keep `--cpu` below the number the number of query organisms.

```
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
                        mycodb, DEFAULT: masterdb
  -o OUTPUT, --output OUTPUT
  -p PREVIOUS, --previous PREVIOUS
                        Previous run directory
  --cpu CPU             DEFAULT: all
```

<br />

## hmmsearch mycoDB
### db2hmmsearch.py
`db2hmmsearch.py` will compile hmmsearch results, optionally ouput fasta/hmmalign to original models/trim alignments from a profile hidden markov model. This script supports multiprocessing and uses all detected cores by default.

e.g.: `db2hmmsearch.py -d profile.hmm` will run the master database referencing the inputted hmmdb. 

`db2hmmsearch.py -d profile.hmm -b 1 -t 75` takes the best hit from hits with 25 - 100% query coverage
```
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

<br />

## Phylogenetic Analysis
### fa2tree.py
`fa2tree.py` will input a fasta file or directory of fasta files, trim, and generate a tree either via job submission (`slurm` or `torque`) or immediate execution. 

e.g. Construct a tree from a fasta
```
fa2tree.py -i <FASTA>.fa -t iqtree
```

<br />

## Hiearchical agglomerative clustering
### aggClus.py
Hiearchical agglomerative clustering is a useful systematic approach to
extracting groups of sequences for phylogenetic analysis. Ubiquitous genes,
like P450s, will often yield 10,000s of results for BLAST searches against the
mycoDB. Constructing and visualizing a tree of this magnitude just is not
practical in many cases; it is therefore necessary to decrease the size of the 
dataset to a workable size while relying on biological information (global
pairwise alignments).

`aggClus.py` will either take a `fasta` and generate a distance matrix using %
identity of `needle` alignments, or optionally work from a tab-delimited
distance file (say from `usearch -calc_distmx`). Subsequently, `aggClus.py`
will cluster sequences via hierarchical agglomerative clustering. Currently,
constructing a distance matrix takes quite long because `needle` cannot
throwout alignments when it detects low % identity in the alignment phase,
whereas `usearch -calc_distmx` does. It is therefore recommended to first
create the distance matrix using `usearch` and then use `aggClus.py` to perform
the agglomerative clustering step instead of `usearch -cluster_aggd`, which
limits dataset size unless the proprietary version is purchased.

e.g. Calculate a distance matrix and cluster from a fasta with a minimum
identity of 0.3 to consider a connection:
```
aggClus.py -f <FASTA>.fa -m 0.7
```

e.g. Construct a usearch distance matrix and use aggClus.py to cluster:
```
usearch -calc_distmx <FASTA>.fa -tabbedout <DISTMX>.tsv -maxdist 0.7
aggClus.py -d <DISTMX>.tsv
```

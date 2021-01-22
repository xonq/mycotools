# mycotools databases
## v0.0.9alpha CURRENTLY NONFUNCTIONAL IF YOU ARE NOT ON OHIO SUPERCOMPUTER

Analyzing large datasets does not have to suck. But how does one acquire the genomes and proteome data; curate to uniformity; efficiently use subsets of the data; update the data; and deal with incorrectly identified accessions in established databases? The mycotools database and associated [scripts](https://gitlab.com/xonq/mycotools_scripts/) are designed for automated curation / maintenance of fungal genomic data from NCBI and JGI and modular downstream analysis of these data. 

mycotools databases are flat, tab delimitted reference files that feauture: proteome, genome, & gene coordinate files; automatic downloads and setup; automated updates; curated gene naming; unique identifiers for each organism (first 3 letters of genus, first 3 of species, unique accession number); taxonomy; publication data; and other relevant metadata. Downstream, my scripts then modularly interface to introduce organisms into an analysis. Via these database tools, one can `abstractDB.py` a database based on the taxonomy, publication information, or organisms relevant for their analysis in one line. Alternatively, one can simply use the default master database to conduct analyses considering all of published Fungi! This modularity is the reason why I created mycotools.

### What is the plan?
I cultivate a master database from all NCBI and JGI genomes through these scripts. I then construct trees from conserved orthologs throughout Fungi and use these data to determine which "fungi" are not actually fungi - a species phylogeny will be uploaded to the repository. From here, I create a list of "forbidden" genome codes and taxonomy IDs that are removed from the master database and will be flagged in this git page. You install my [scripts](https://gitlab.com/xonq/mycotools_scripts/) and run `updateDB.py --init` to initialize the database install based on my reference database. Then get going! All you have to do to update the database is simply run `updateDB.py -y` to update whenever available. The goal is to recreate the trees each year and have a "stable" branch that only includes organisms evaluated through the trees, while also having an "unstable" branch that includes tree-evaluated organisms and the most recent accessions that have not been incorporated into a tree.

<br /><br />

# Use
Because these databases are flat, tab delimitted files, they can be queried and manipulated using bash / Unix commands or alternatively loaded into Excel - provided the output is saved as a plain text `.tsv` without hidden characters.

<br />

# Scripts
`abstractDB.py`: abstracts a sub database based on taxonomy, publication information, source, or genome codes

`jgi2db.py`: downloads JGI genomes and assimilates into a database or prepares an update for an existing database

`ncbi2db.py`: downloads NCBI genomes and assimilates into a database or prepares an update for an existing database

`updateDB.py`: updates database using a prepared reference database, or simply updates the existing information

<br /><br />

# Install
#### Preparing
Make sure you have a working python3 environment, then install the dependencies:
```
pip install mycotools Biopython pandas
```
Next, use a text editor to add where you will install database information to your `.bash_profile` or `.bashrc`:
```
export DB=<DB/INSTALL/PATH>
export PROTEOME=$DB/proteome
export ASSEMBLY=$DB/assembly
export GFF3=$DB/gff3
export BLAST=$PROTEOME/blastdb
```

<br />

#### Set up
I have plans to create a `setupDB.py` script in the future, but for now, clone this repository and use `jgi2db.py` followed by `ncbi2db.py` then `updateDB.py` to get started making your own mycotools database. These scripts will create a database file, download genomic data, curate the data, acquire metadata, and assemble them into `$DB`.

<br />

#### Updating
Updating is easy: simply rerun the `ncbi2db.py` and `jgi2db.py` scripts followed by `updateDB.py`. JGI genomes by default take precedence over detected duplicates in NCBI.

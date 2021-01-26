# PURPOSE
MycoTools is a compilation of computational biology tools and database (MycoDB) software designed to increase throughput analyzing fungal genomic data (JGI & NCBI). MycoDB is a database schema with uniform gene coordinate (`gff`) and fasta header curation, systematic code naming, and modularity as guiding principles, which enable analyses with datasets of 1000s of fungi to as few as 1 or 2. 

By integrating with the curated MycoDB, MycoTools cuts time on routine tasks like retrieving `gff` or `fasta` accessions, running and compiling `fasta`s of MycoDB BLAST/hmmsearches, and other analyses like creating a phylogeny with `fa2tree.py` or wrangling your dataset with agglomerative clustering via `aggClus.py`. MycoTools includes sets of utilities that also enable easy acquisition of batches of sequence data using `ncbiDwnld.py` and `jgiDwnld.py`. Check out the [USAGE.md](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md) for a guide on the possibilities MycoTools can enable in your research.

<br />

If MycoTools and/or MycoDB contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the MycoTools version in line.

---

<br />

# UPDATE
```
pip install mycotools --upgrade
```

<br />

# INSTALL
## 1) Install miniconda3
Miniconda3 is an environment manager that allows you to install your own `python` and isolate it from the system `python`. Type `yes` at the end to run `conda init`. Pay attention to the installation if you want to install to a specific path (e.g. `~/software/miniconda3`- make sure to include `miniconda3`). This keeps your home folder from getting cluttered.

```	
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > ~/miniconda3.sh
bash ~/miniconda3.sh
```

If you are using the Ohio Supercomputer, or are having problems with initializing `conda`, run this to autostart your conda environment upon login:
```
echo "\nsource ~/.bashrc" >> ~/.bash_profile
```

<br />

## 2) Install MycoTools
Restart your login. If upon login you do not see `(base)` to the left of your bash shell prompt then *anytime you want to use the python environment you created, you have to activate it via* `source activate base`. If necessary, activate the environment, then use `pip` to install MycoTools:

```
pip install mycotools
```

<br />

## 3) Install MycoDB 
#### OSC
Many of my scripts interface with [MycoDBs](https://gitlab.com/xonq/mycodb/-/blob/master/README.md). If you are using the Ohio Supercomputer and have access to PAS1046, then simply append these commands to your `~/.bash_profile` in your OSC home folder.
```
export MYCODB=/fs/project/PAS1046/databases/konkel 	# database
export MYCOFNA=$MYCODB/assembly/ 	# database assemblies
export MYCOFAA=$MYCODB/proteome/ 	# database proteomes
export MYCOGFF3=$MYCODB/gff3/ 		# database gff3s
```

Restart your login and you're good to proceed to the [usage guide!](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md)

#### Non-OSC
If you are not on OSC or want to learn more about MycoDB, you can copy your own MycoDB by following [this guide](https://gitlab.com/xonq/mycodb/-/blob/master/README.md). Essentially run `updateDB.py --init <NEW_DATABASE_DIRECTORY>` to get started with default options.

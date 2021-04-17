# PURPOSE
MycoTools is a compilation of computational biology tools and database (MycoDB) software designed to increase throughput analyzing fungal genomic data (JGI & NCBI). MycoDB is a database schema with uniform gene coordinate (`gff`) and fasta header curation, systematic code naming, and modularity as guiding principles, which enable analyses with datasets of 1000s of fungi to as few as 1 or 2. MycoDB is not available outside of the Ohio SuperComputer until we publish the tools. Please email `konkelzach@protonmail.com` if you are interested in becoming an early adopter.

By integrating with the curated MycoDB, MycoTools cuts time on routine tasks like retrieving `gff` or `fasta` accessions, running and compiling `fasta`s of MycoDB BLAST/hmmsearches, and other analyses like creating a phylogeny with `fa2tree.py` or wrangling your dataset with agglomerative clustering via `aggClus.py`. MycoTools includes sets of utilities that also enable easy acquisition of batches of sequence data using `ncbiDwnld.py` and `jgiDwnld.py`. Check out the [USAGE.md](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md) for a guide on the possibilities MycoTools can enable in your research. 

<br />

If MycoTools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the MycoTools version in line.

---

<br />

# UPDATE
```
pip install mycotools --upgrade
```

<br />

# INSTALL
## 1) Install miniconda3
Miniconda3 is an environment manager that allows you to install your own `python` and isolate it from the system `python`. If you are using the Ohio Supercomputer Center (OSC) you can access conda by creating an environment as described at the top of [this link](https://www.osc.edu/resources/getting_started/howto/howto_add_python_packages_using_the_conda_package_manager). Then proceed to step 2.
Otherwise, pay attention to the installation if you want to install to a specific path (e.g. `~/software/miniconda3`- make sure to include `miniconda3`). This keeps your home folder from getting cluttered. 

```	
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > ~/miniconda3.sh
bash ~/miniconda3.sh
```

Restart your shell and create a new conda environment
```
conda create -n mycotools python=3.8
conda activate mycotools
```

<br />

## 2) Install MycoTools
With your environment activated
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
MycoDB is currently not available for widespread use. We are waiting to complete our intial analyses and publish the data before releasing.

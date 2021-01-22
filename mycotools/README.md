# PURPOSE
These scripts are a compilation of computational biology tools to increase analysis throughput. In particular, these scripts are tailored toward use with the `mycotools database`, a database schema with modularity as the guiding principle, enabling analyses with datasets of 1000s of fungi to as few as 1 or 2.

---

<br />

# INSTALLING SCRIPTS
## 1) Download & install miniconda3
Miniconda3 is an environment manager that allows you to install your own `python` and isolate it from the system `python`. Pay attention to the installation if you want to install to a specific path (e.g. `~/software/miniconda3`- make sure to include `miniconda3`). This keeps your home folder from getting cluttered.

```	
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > ~/miniconda3.sh
bash ~/miniconda3.sh
```

If you are using the Ohio Supercomputer, or are having problems with initializing `conda`, edit and add this line to your `~/.bash_profile`:
```
export PATH="<ABSOLUTE/PATH/MINICONDA/INSTALL>/bin:$PATH"
```

<br />

## 2) Install scripts & dependencies
Restart your login. *Anytime you want to use the python environment you just installed, you may have to activate it via* `source activate base`. There are relatively easy ways to [make this automatic](https://docs.anaconda.com/anaconda/user-guide/faq/), but it is beyond the scope of this install. 

```
source activate base
pip install mycotools Biopython pandas scikit-learn
```

<br />

## 3) Installing mycodb 
#### OSC
Many of my scripts interface with [mycotools databases](https://gitlab.com/xonq/scripts/-/blob/master/database/README.md). If you are using the Ohio Supercomputer and have access to PAS1046, then simply append these commands to your `~/.bash_profile` in your OSC home folder.
```
export MYCODB=/fs/project/PAS1046/databases/konkel 	# database
export MYCOFNA=$MYCODB/assembly/ 	# database assemblies
export MYCOFAA=$MYCODB/proteome/ 	# database proteomes
export MYCOGFF=$MYCODB/gff/ 		# database gffs
export MYCOGFF3=$MYCODB/gff3/ 		# database gff3s
export MYCOBLAST=$MYCODB/proteome/blastdb 	# blastdbs
export MYCOSMASH=/fs/project/PAS1046/projects/secmet/antismash 	# antismash outputs
```

Restart your login and you're good!

#### Non-OSC
If you are not on OSC or want to learn more about mycotools, you can create your own mycotools database following [this guide](https://gitlab.com/xonq/mycotools_scripts/-/blob/master/database/README.md).

These scripts are currently under heavy development. Before running analyses, it is a good idea to check for upgrades by running `pip install mycotools --upgrade` after activating your conda environment.

These are scripts designed to automate many processes in my computational biology work and interface with the mycotools database. If you want the best experience with them, please install the database as described after the script installation. 

<br />

# INSTALLING SCRIPTS
## 1) Download & install miniconda3
You will need a virtual python3 environment, which essentially is installing you your own local python and isolating it in its own "environment" from the python your system relies on to operate - this allows you to install python packages without interfering with the system and is required to install packages on HPCs like Ohio Supercomputer (OSC). All you need to do is pay attention to the installation if you want to install to a specific path (I install mine to a software folder in my home directory, e.g. `~/software/miniconda3` - make sure to add `miniconda3` as the final portion of the path). Otherwise, answer `y` to all:

```	
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > ~/miniconda3.sh
bash ~/miniconda3.sh
```

Then edit and add this line to your `~/.bash_profile`:
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

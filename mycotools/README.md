# PURPOSE
MycoTools is a compilation of computational biology tools and database (MycoDB) maintenance software directed toward the analysis of fungal sequence data. MycoDB is a database schema with modularity as the guiding principle, enabling analyses with datasets of 1000s of fungi to as few as 1 or 2.

<br />

# CITING
If you use MycoTools or the MycoDB in your analyses, please cite this git repository (gitlab.com/xonq/mycotools) and mention MycoTools in line until I have published these scripts in the literature.

---

<br />

# INSTALL
## 1) Install miniconda3
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

## 2) Install MycoTools
Restart your login. If upon login you do not see `(base)` to the left of your bash shell prompt then *anytime you want to use the python environment you created, you have to activate it via* `source activate base`. There are relatively easy ways to [make this automatic](https://docs.anaconda.com/anaconda/user-guide/faq/), but it is beyond the scope of this install. If necessary, activate the environment, then use `pip` to install MycoTools:

```
pip install mycotools
```

<br />

## 3) Install MycoDB 
#### OSC
Many of my scripts interface with [mycotools databases](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/db/README.md). If you are using the Ohio Supercomputer and have access to PAS1046, then simply append these commands to your `~/.bash_profile` in your OSC home folder.
```
export MYCODB=/fs/project/PAS1046/databases/konkel 	# database
export MYCOFNA=$MYCODB/assembly/ 	# database assemblies
export MYCOFAA=$MYCODB/proteome/ 	# database proteomes
export MYCOGFF3=$MYCODB/gff3/ 		# database gff3s
```

Restart your login and you're good to proceed to the [usage guide!](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md)

#### Non-OSC
If you are not on OSC or want to learn more about mycotools, you can create your own mycotools database following [this guide](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/db/README.md). Essentially run `updateDB.py --init <NEW_DATABASE_DIRECTORY>` to get started with default options.

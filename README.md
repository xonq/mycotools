<p align="center">
    <img
        src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/pictogo.white.png"
    >
</p>

<br /><br />

# NOTE
This software is a beta release (prokaryote is alpha state) - errors are expected. Kindly report them.
If you can find the bug, even better! The goal is to reach a longterm stable
release, though maintaining the software for my use is currently the priority.

# PURPOSE
Bring broadscale comparative genomics to the masses. 

Mycotools is a compilation of computational biology tools and database [MycotoolsDB](https://github.com/xonq/mycotools/blob/master/MTDB.md) software that facilitate large-scale fungal comparative genomics. MycotoolsDB locally assimilates all NCBI and MycoCosm (Joint Genome Institute) genomes into a database schema with uniform file curation, scalability, and automation as guiding principles. 

- Database initialization is as simple as `mtdb u --init <DIR>`
- `mtdb u --update` brings the database to the current date
- The MycotoolsDB (MTDB) uniformly curates the numerous iterations of
  the `gff`, allowing for reliable analyses and format expectations from
  multiple eras
- The `.mtdb` database format enables swift transitions from analyses with datasets of 100,000s genomes to as few as a lineage of interest
- The MycotoolsDB can be adminstered by one administrator and accessed by
  multiple users that integrate via `mtdb -i <DB_DIR>`

<p align="center">
    <img
        src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/mtdb.png"
    >
</p>

<br />

Mycotools facilitates routine-complex
tasks like retrieving locus, `gff`, or `fasta` accessions; running and compiling
`fasta`s of MycotoolsDB BLAST/hmmsearches; automated phylogenetic analysis
pipelines from BLAST to Pfam extraction to tree prediction, etc etc. Mycotools
includes sets of utilities that also enable easy acquisition of batches of
sequence data using `ncbiDwnld.py` and `jgiDwnld.py`. Please see the [USAGE
guide](https://github.com/xonq/mycotools/blob/master/USAGE.md) for
more information.

<br />

# USAGE
Check out [README.md](https://github.com/xonq/mycotools/blob/master/README.md) for install and the [USAGE.md](https://github.com/xonq/mycotools/blob/master/USAGE.md) for a guide. 

<br />

# CITING

If Mycotools contribute to your analysis, please cite this git repository (github.com/xonq/mycotools) and mention the Mycotools version in line.

---

<br />

# INSTALL

## 1. Installing miniconda
Miniconda3 is a software environment manager:

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > ~/miniconda3.sh
bash ~/miniconda3.sh
```

Activate miniconda and initialize it so it starts up automatically
```bash
source activate base # if this fails run conda activate base
conda init
```


<br />

## 2. Setting up miniconda
Setup and prioritize channels for your miniconda installation. This step must be
completed for new and old installs.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

<br />

## 3. Installing mycotools
Make sure `conda` is active, usually by seeing `(base)` in in your shell.
If not, try `conda activate base` or `source activate base`. 

```bash
conda create -n mycotools mycotools -c xonq -y
conda activate mycotools
python3 -m pip install mycotools --upgrade
mtdb -d
```

NOTE: There is an unfortunate build dependency conflict with the conda package,
so be sure to explicitly run `pip install` as noted above to update the
mycotools package inside the conda environment.

<br />

Determine if you are going to link to an already installed database, or become
the administrator of a new one:

## 4a. USER: Integrate with already initialized MycotoolsDB
To link with an existing database, fill in `<PATH>` with the database path

```bash
mtdb --init <DB_PATH>
```

You're good to proceed to the
[usage guide!](https://gitlab.com/xonq/mycotools/-/blob/master/USAGE.md)

## 4b. ADMINISTRATOR: Initialize a local MycotoolsDB
```bash
mtdb update --init <DB_PATH>
```

<br />

# UPDATE
Mycotools is currently in an advanced beta state with frequent updates. It is
recommended to run the following in your conda environment if you are having
trouble with analyses:

```bash
python3 -m pip install mycotools --upgrade
```

NOTE: Make sure the conda environment is active when updating.
I recommend updating with `pip` because the `conda` distribution 
is currently prone to dependency issues and will not update reliably.


<br /><br /><br />

### A NOTE ON THE CODE
Each standalone script is written with `__name__ == '__main__'`, designed to
handle running the script from the command line, as well as `main` function(s),
which are importable modules executing the purpose of the script. This enables Mycotools
to be a pipelining-friendly software suite, both from a command line and
python scripting standpoint.

Code edits should focus on stabilizing existing features and simplifying/decerasing the code base.
I try to implement code aligned with principles of the [functional
programming paradigm](https://docs.python.org/3/howto/functional.html) and
modifications should act in accord with this paradigm, i.e. sparing
implementation of new classes, limited necessary abstraction, no hidden state
changes, and function-based flow.


<img align="right" src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/ablogo.png">

<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />

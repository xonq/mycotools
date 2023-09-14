<p align="center">
    <img
        src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/pictogo.white.png"
    >
</p>

<br /><br />

# NOTE
This software is a beta release (prokaryote is alpha state) - errors are
expected. Kindly raise an issue if the error has not been mentioned.
If you can find the bug, even better! The goal is to reach a longterm stable
release with well-annotated PEP8 code, though maintaining the software for my 
use is currently the priority.

# PURPOSE
Bring broadscale comparative genomics to the masses. 

Mycotools is a compilation of computational biology tools and database
([MycotoolsDB](https://github.com/xonq/mycotools/blob/master/MTDB.md)) software that facilitate large-scale comparative genomics. MycotoolsDB locally assimilates all NCBI and MycoCosm (Joint Genome Institute) genomes into a database schema with uniform file curation, scalability, and automation as guiding principles. 

- Database initialization is as simple as `mtdb u --init <DIR>`
- `mtdb u --update` brings the database to the current date
- The MycotoolsDB (MTDB) uniformly curates the numerous iterations of
  the `gff` into a consistent standard
- The `.mtdb` database format enables swift transitions from analyses with datasets of 100,000s genomes to as few as a lineage of interest
- Mycotools software suite includes modules to automate routine-complex
  comparative genomics

<p align="center">
    <img
        src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/mtdb.png"
    >
</p>

<br />

# CITING

If Mycotools contribute to your analysis, please cite the preprint and mention
the version in-line. 

Konkel, Z., Slot, J. C. Mycotools: An Automated and Scalable Platform for
Comparative Genomics. bioRxiv 2023.09.08.556886; doi: https://doi.org/10.1101/2023.09.08.556886

---

<br />

# INSTALL

The installation guide will use miniconda3 as the environment manager. Please
reference their [install and initialization instructions](https://docs.conda.io/projects/miniconda/en/latest/).

<br />

## 1. Configuring miniconda3
Setup and prioritize channels for your miniconda installation. This step must be
completed for new and old installs.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

<br />

## 2. Installing mycotools
Make sure `conda` is active, usually by seeing `(base)` in in your shell.
If not, try `conda activate base` or `source activate base`. 

```bash
conda create -n mycotools mycotools -c xonq -y
conda activate mycotools
python3 -m pip install mycotools --upgrade
mtdb -d
```

NOTE: The conda package is broken on its own, so please 
explicitly run `pip install` as noted above to update the
mycotools package inside the conda environment.

<br />

Determine if you are going to link to an already installed database, or become
the administrator of a new one:

## 3a. USER: Integrate with already initialized MycotoolsDB
To link with an existing database, fill in `<PATH>` with the database path

```bash
mtdb --init <DB_PATH>
```

<br />

## 3b. ADMINISTRATOR: Initialize a local MycotoolsDB
```bash
mtdb update --init <DB_PATH>
```

<br /><br />

# USAGE

Once installed, you're good to proceed to the
[usage guide!](https://github.com/xonq/mycotools/blob/master/USAGE.md)


<br /><br />

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
Each standalone script is written with a `cli` function, designed to
handle running the script from the command line, as well as `main` function(s),
which are importable python modules. This enables Mycotools
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

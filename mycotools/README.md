![mycotools](https://gitlab.com/xonq/mycotools/-/raw/master/misc/pictogo.png)

<br /><br />

# CITING

If Mycotools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the Mycotools version in line.

---

<br />

# UPDATING
Mycotools is currently in an advanced beta state with frequent updates. It is
recommended to run the following in your conda environment before analyses:

```bash
conda update mycotools -c xonq
```

<br />

# INSTALL
## Installing miniconda
Miniconda3 is a software environment manager. Pay attention to the installation if you want to install to a specific path (e.g. `~/software/miniconda3`).

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

## Setting up miniconda
Setup and prioritize channels for your miniconda installation. This step must be
completed for new and old installs.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

<br />

## Installing mycotools
Create a new conda environment; `(base)` should appear
in your shell before running - if not, try `conda activate base` or `source activate base`.

```bash
conda create -n mycotools mycotools -c xonq
conda activate mycotools
mtdb -d
```

NOTE: if you are having trouble on install, run the `conda config` commands
above

<br />

## Integrate with installed MycotoolsDB 
#### OSC
If you are using the Ohio Supercomputer and have access to PAS1046/PAS1568, you can use the preinstalled fungal and prokaryote databases respectively without reinitializing. Otherwise, change the `<PATHS>` to your database administrators path. Initializing databases is covered in the [usage guide](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md).
Run this command with the respective path to integrate your installation with the fungal or prokaryote database:
```bash
mtdb --init <PATH>
```
Fungi: `/fs/project/PAS1046/databases/mycotoolsdb/`
Prokaryote: `/fs/ess/PAS1568/mycotools/`


You're good to proceed to the [usage guide!](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md)

#### Non-OSC
MycotoolsDB will be available September, 2022.

<br /><br /><br />

### A NOTE ON THE CODE
Each standalone script is written with `__name__ == '__main__'`, designed to
handle running the script from the command line, as well as `main` function(s),
which are importable modules executing the purpose of the script. This enables Mycotools
to be a pipelining-friendly software suite, both from a command line and
python scripting standpoint.

I *primarily* abide by the [functional
programming paradigm](https://docs.python.org/3/howto/functional.html).
Furthermore, I only create new class objects if a task *needs specific manipulations* that are cumbersome in default classes alone. Any code edits should follow this guideline and adhere to the functional paradigm i.e. no hidden state changes.


<img align="right" src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/ablogo.png">

<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />

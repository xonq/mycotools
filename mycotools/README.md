![mycotools](https://gitlab.com/xonq/mycotools/-/raw/master/misc/pictogo.png)

<br /><br />

# CITING

If Mycotools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the Mycotools version in line.

---

<br />

# UPDATE
Mycotools is currently in an advanced beta state with frequent updates. It is
recommended to run the following in your conda environment if you are having
trouble with analyses:

```bash
conda update mycotools -c xonq
```

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
mtdb -d
```

<br />

## 4. Integrate with installed MycotoolsDB 
If you are using the Ohio Supercomputer and have access to PAS1046/PAS1568, you can use the preinstalled fungal and prokaryote databases. If you are outside OSC change the `<PATHS>` to your database administrators path. Initializing databases is covered in the [usage guide](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md).

To link with an existing database, fill in `<PATH>` with the database path. OSC databases denoted below.
```bash
mtdb --init <PATH>
```
Fungi: `/fs/project/PAS1046/databases/mycotoolsdb/`
Prokaryote: `/fs/ess/PAS1568/mycotools/`


You're good to proceed to the [usage guide!](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md)

#### Non-OSC
MycotoolsDB will be available May, 2023.

<br /><br /><br />

### A NOTE ON THE CODE
Each standalone script is written with `__name__ == '__main__'`, designed to
handle running the script from the command line, as well as `main` function(s),
which are importable modules executing the purpose of the script. This enables Mycotools
to be a pipelining-friendly software suite, both from a command line and
python scripting standpoint.

Code edits should focus on stabilizing existing features and simplifying/decerasing the code base.
I try to implement code aligned with principles of the [functional
programming paradigm](https://docs.python.org/3/howto/functional.html). 


<img align="right" src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/ablogo.png">

<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />

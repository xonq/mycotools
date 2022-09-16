![mycotools](https://gitlab.com/xonq/mycotools/-/raw/master/misc/pictogo.png)

<br /><br />

# CITING

If Mycotools contribute to your analysis, please cite this git repository (gitlab.com/xonq/mycotools) and mention the Mycotools version in line.

---

<br />

# UPDATE
```
pip install mycotools --upgrade
```

<br />

# INSTALL
## Installing miniconda and bioconda
Miniconda3 is an environment manager that allows you to install your own `python` and isolate it from the system `python`. If you are using the Ohio Supercomputer Center (OSC) you can access conda by creating an environment as described at the top of [this link](https://www.osc.edu/resources/getting_started/howto/howto_add_python_packages_using_the_conda_package_manager). Then proceed to step 2.
Otherwise, pay attention to the installation if you want to install to a specific path (e.g. `~/software/miniconda3`- make sure to include `miniconda3`). This keeps your home folder from getting cluttered. 

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > ~/miniconda3.sh
bash ~/miniconda3.sh
```

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

<br />

## Installing mycotools
Restart your shell and create a new conda environment; `(base)` should appear
in your shell (if not, try `conda activate base` or `source activate base`):
```bash
conda create -n mycotools mycotools -c xonq
conda activate mycotools
```

NOTE: if you are having trouble on install, run the `conda config` commands
mentioned above; if you had to run `conda activate base`/`source activate base`
then you will need to add the respective command to your profile, e.g. `echo
"<CMD> activate base" >> ~/.bash_profile`, to automatically call conda.

<br />

## Integrate with installed MycotoolsDB 
#### OSC
If you are using the Ohio Supercomputer and have access to PAS1046/PAS1568, then simply run this command with the respective path to integrate your installation with the fungal or prokaryote database:
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
python scripting standpoint,  while also emphasizing the Unix
philosophy for each script to *do one thing and do it well*. 

I *primarily* abide by the [functional
programming paradigm](https://docs.python.org/3/howto/functional.html).
I only create new class objects if 1) used frequently, 2) requires rigid formatting, 
and 3) *needs specific manipulations* that are cumbersome in default classes alone. 
Any code edits should follow this guideline and implement the functional paradigm.


<img align="right" src="https://gitlab.com/xonq/mycotools/-/raw/master/misc/ablogo.png">

<br /><br /><br /><br /><br /><br /><br /><br /><br /><br /><br />

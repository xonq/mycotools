![mycotools](https://gitlab.com/xonq/mycotools/-/raw/master/misc/logo.png)

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
## 1) Install miniconda3
Miniconda3 is an environment manager that allows you to install your own `python` and isolate it from the system `python`. If you are using the Ohio Supercomputer Center (OSC) you can access conda by creating an environment as described at the top of [this link](https://www.osc.edu/resources/getting_started/howto/howto_add_python_packages_using_the_conda_package_manager). Then proceed to step 2.
Otherwise, pay attention to the installation if you want to install to a specific path (e.g. `~/software/miniconda3`- make sure to include `miniconda3`). This keeps your home folder from getting cluttered. 

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > ~/miniconda3.sh
bash ~/miniconda3.sh
```

Restart your shell and create a new conda environment; `(base)` should appear
in your shell:
```bash
conda create -n mycotools python=3.8
conda activate mycotools
```

<br />

## 2) Install Mycotools
With your environment activated (`(mycotools)` should appear in your shell):
```bash
pip install mycotools
```

To upgrade to the latest version:
```bash
pip install mycotools --upgrade
```

<br />

## 3) Integrate with installed MycotoolsDB 
#### OSC
If you are using the Ohio Supercomputer and have access to PAS1046, then simply run this command to integrate your installation with the eukaryote database:
```bash
mycodb --init /fs/project/PAS1046/databases/mycodb/
```

for the prokaryote OSC database:
```bash
mycodb --init /fs/ess/PAS1568/mycotools
```

Restart your login and you're good to proceed to the [usage guide!](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md)

#### Non-OSC
MycotoolsDB is currently not available for widespread use. We are waiting to complete our intial analyses and publish the data before releasing.

<br /><br /><br />

### A NOTE ON THE CODE
Each standalone script is written with `__name__ == '__main__'`, designed to
handle running the script from the command line, as well as `main` function(s),
which are importable modules that run the analysis. This enables Mycotools
to be a pipelining-friendly software suite, both from a command line and
python scripting standpoint,  while also emphasizing the Unix
philosophy for each script to *do one thing and do it well*. 

I have found that I am *primarily* a [functional
programmer](https://docs.python.org/3/howto/functional.html), i.e. my code seldom incorporates novel classes.
I am not a purist, as you will find the `mtdb` class for the 
Mycotools Database objects (in the future there will be several more classes),
but only because it is 1) used frequently, 2) requires rigid formatting, and 3) *needs
specific manipulations* that are cumbersome in default classes alone. 
*I will continue functional programming*, striving to format 
in accord with PEP-8 with minimal deviance - while there is still work to be done any code edits should 
favor the functional paradigm, unless a demonstrated need for a class (above) exists.

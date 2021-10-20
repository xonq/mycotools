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

## 2) Install Mycotools
With your environment activated
```
pip install mycotools
```

<br />

## 3) Install MycotoolsDB 
#### OSC
Many of my scripts interface with [MycotoolsDBs](https://gitlab.com/xonq/MycotoolsDB/-/blob/master/README.md). If you are using the Ohio Supercomputer and have access to PAS1046, then simply run this command to integrate your installation with the database:
```
mycodb --init /fs/project/PAS1046/databases/mycodb/
```

Restart your login and you're good to proceed to the [usage guide!](https://gitlab.com/xonq/mycotools/-/blob/master/mycotools/USAGE.md)

#### Non-OSC
MycotoolsDB is currently not available for widespread use. We are waiting to complete our intial analyses and publish the data before releasing.

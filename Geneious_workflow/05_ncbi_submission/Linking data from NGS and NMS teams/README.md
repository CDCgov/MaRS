## README

Please set up a virtual env for all python jupyter notebooks. Read this article [Use VE insdie Jupyter Notebook and Lab - Best Practices](https://www.zainrizvi.io/blog/jupyter-notebooks-best-practices-use-virtual-environments/)

The [venv](https://docs.python.org/3/library/venv.html) module allows for creating lightweight virtual environments with their own site directories, optionally siolated from system directories. Each VE has its own Python binary and independent set of python packages.

If you have python 3.x installed, [venv](https://docs.python.org/3/library/venv.html) is already installed as well.

### Creating a virtual env (ve)

Start up `jupyter lab` in current working directory and open up a terminal. The below commands will set up a _virtual environment_ and show up in the jupyter launcher.  


```python
## If in current desired directory ##

python3 -m venv data_link_ve --system-site-packages      # create data_link_ve env

source data_link_ve/bin/activate  # activate the ve
```
If correctly initialized, `(data_link_ve)` should appear now in front of the `userID$`

### Registering this ve with jupyter lab

```python
## It should show up in launcher list once page is refreshed ##

python -m ipykernel install --user --name=data_link_ve
```

### To install specific packages or modules

Use the terminal to install any specific modules and dependencies in this _virtual environment_. However, if for some reason you install them system wide, they should still be able to be called because when we specified `--system-site-packages` during intial creation of the _virtual environment_. 

```python3
python3 -m pip install `module name`

```
**Any python modules or packages installed now via pip will only persist in this virtual environment**

### Created virtual env directories

Once initialized, the _virtual environment_ will create a new directory with all installed dependencies based on the specified `$PATH` when the _ve_ was created. Using the above example, a directory called `data_link_ve` should have been created.

### Exiting the virtual environment

```python
## Make sure to deactivate the ve once done working in jupyter ##
deactivate  # exit the virtual env

```

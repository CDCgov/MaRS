> Author: @ET 4/5/22 :goat:
>> Edited: N/A

### Backgorund ###

 Read this article [Use VE insdie Jupyter Notebook and Lab - Best Practices](https://www.zainrizvi.io/blog/jupyter-notebooks-best-practices-use-virtual-environments/). A `virtual environment` will need to be set up for the python jupyter notebooks. The [venv](https://docs.python.org/3/library/venv.html) module allows for creating lightweight virtual environments with their own site directories, optionally siolated from system directories. Each VE has its own Python binary and independent set of python packages.

### Python set up ###
If you have python 3.x installed, [venv](https://docs.python.org/3/library/venv.html) is already installed as well. Open up a terminal and run `which python`. This will print the python version you have installed. If you get nothing, no worries. We need to add python version management using [pyenv](https://github.com/pyenv/pyenv). For Unix OS, we also need install the package manager [homebrew](https://brew.sh/).

```bash
# Install home brew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install pyenv and add to .zshrc
brew install pyenv.      
brew install 3.10.1        # change to latest version; this was 3.10.1 on 4/5/22   
pyenv global 3.10.1    
pyenv version              # check it worked and correct version is shown

# add to shell profile to load up on start (change to .bash_profile on Linux)
echo -e $'if command -v pyenv 1>/dev/null 2>&1; then\n export PYENV_ROOT="$HOME/.pyenv"\n export PATH="$PYENV_ROOT/bin:$PATH"\n eval "$(pyenv init --path)"\n eval "$(pyenv init -)"\nfi' >> ~/.zshrc
```

### Creating a virtual env (ve)

Start up `jupyter lab` in current working directory and open up a terminal. The below commands will set up a _virtual environment_ and show up in the jupyter launcher.  


```bash
## If in current desired directory ##  

python3 -m venv jupyter_ve --system-site-packages      # create jupyter_ve

source jupyter_ve/bin/activate                         # activate the ve
```
If correctly initialized, `(jupyter_ve)` should appear now in front of the `userID$`

### Registering this ve with jupyter lab

```bash
## It should show up in launcher list once page is refreshed ##

python -m ipykernel install --user --name=jupyter_ve
```

### To install specific packages or modules

Use the terminal to install any specific modules and dependencies in this _virtual environment_. However, if for some reason you install them system wide, they should still be able to be called because when we specified `--system-site-packages` during intial creation of the _virtual environment_.

```bash
python3 -m pip install <module name>

```
**Any python modules or packages installed now via pip will only persist in this virtual environment**

### Created virtual env directories

Once initialized, the _virtual environment_ will create a new directory with all installed dependencies based on the specified `$PATH` when the _ve_ was created. Using the above example, a directory called `jupyter_ve` should have been created.

### Exiting the virtual environment

```bash
## Make sure to deactivate the ve once done working in jupyter ##
deactivate    # exit the virtual env

```

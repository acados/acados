# acados documentation

Based on `sphinx` and markdown

### Prerequisites

Get and install doxygen:
```
wget http://doxygen.nl/files/doxygen-1.8.15.linux.bin.tar.gz
tar xf doxygen-1.8.15.linux.bin.tar.gz
# add path to doxygen/bin to you bashrc
```

This guide shows the *easy* way of using conda to manage virtual environments. You can of course use the environment of your choice to install pip packages.

Get and install miniconda (follow the installation instructions, best to keep default path ~/miniconda3, and say yes to initialize miniconda in your bashrc):
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc # or restart your terminal
```

Create an environment for python and activate the environment
```
conda create -n acados_doc python=3.7
conda activate acados_doc
```

Create a python virtualenv and install `sphinx` and `recommonmark`

```
conda install sphinx=1.8 # or pip install sphinx
pip install recommonmark
pip install breathe
pip install sphinx_rtd_theme
```

### Building the webpage

After doing your modification to the webpage, you can build the webpage
in the folder with the following command:

```
# For Linux, MacOS
cd docs
make html
```

```
# For Windows
cd docs
make.bat html
```

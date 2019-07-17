# acados documentation

Based on `sphinx` and markdown

## Linux

### Prerequisites

* Get and install doxygen and graphviz:
```
sudo apt-get install doxygen graphviz
```
### Python virtual environment

* Get python3.5 or later if not already present.

<!-- ####  -->
<!-- Get and install miniconda (follow the installation instructions, best to keep default path ~/miniconda3, and say yes to initialize miniconda in your bashrc): -->
<!-- ``` -->
<!-- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -->
<!-- sh Miniconda3-latest-Linux-x86_64.sh -->
<!-- source ~/.bashrc # or restart your terminal -->
<!-- ``` -->

<!-- Create an environment for python and activate the environment -->
<!-- ``` -->
<!-- conda create -n acados_doc python=3.7 -->
<!-- conda activate acados_doc -->
<!-- ``` -->
<!-- conda install sphinx=1.8 # or pip install sphinx -->

<!-- Create a python virtualenv and install `sphinx` and `recommonmark` -->

<!-- ``` -->
<!-- pip install -r requirements.txt -->
<!-- ``` -->

* Optionally create and activate a python virtual environment

```
python3 -m venv /path/to/new/virtual/environment
source /path/to/new/virtual/environment/bin/activate
```

* Install requirements:
```
pip3 install -r requirements.txt
```

### Building the webpage

After doing your modification to the webpage, you can build the webpage
in the folder with the following command:

```
# For Linux, MacOS
cd docs
make
```

```
# For Windows
cd docs
make.bat
```

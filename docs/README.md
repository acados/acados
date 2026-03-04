# acados documentation

Based on `sphinx` and markdown


## Prerequisites

### Linux

* Get and install doxygen and graphviz:
```
sudo apt-get install doxygen graphviz
```

### Python virtual environment

* Get python3.8 or later if not already present.
* NOTE: the latest usable version of python is 3.11 as several library modules used by `sphinx 6.0.0.` were deprecated and subsequently removed (eg. `imghdr`).
* Optionally create and activate a python virtual environment (LINUX):
```
python3 -m venv <path_to_virtual_env> # you can use path_to_virtual_env = "env"
source <path_to_virtual_env>/bin/activate
```
* Optionally create and activate a python virtual environment (WINDOWS):
```
python -m venv venv
.\venv\Scripts\Activate.ps1
```

NOTE for Windows: PowerShell scripts are not allowed by default, using [this guide](https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_execution_policies?view=powershell-5.1) you can choose to allow them permanently or only for one session.

* Install requirements:
```
pip install -r requirements.txt
pip install -e <acados_root>/interfaces/acados_template
```

## Build the webpage

To build the website run the following command:

```
# For Linux, MacOS
cd docs
make
```

```
# For Windows
cd docs
.\make.bat html
```

NOTE: The `_build_doxygen_c_interface` folder can have a different auto-generated name and path to the `index.xml` file. This then has to be modified in the `conf.py` configuration file located in `acados\docs` in line 75. (example of the path on Windows11: `\acados\docs\doxygen\_build_doxygen\xml` with the required configuration in line 75 as: `breathe_projects = { "acados": "doxygen/_build_doxygen/xml/" }`).

## Check the generated page
To check the generated html run:
```
# For Linux, MacOS
firefox _build/index.html
```

```
# For Windows - using the default browser
start _build/index.html
```


## Upload the generated files to the server
To upload the webpage you built to syscop.de, run:
```
# For Linux, MacOS -- your ssh key must be authorized on the server.
make upload
```

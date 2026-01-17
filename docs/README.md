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
* Optionally create and activate a python virtual environment
```
python3 -m venv <path_to_virtual_env> # you can use path_to_virtual_env = "env"
source <path_to_virtual_env>/bin/activate
```

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
make.bat
```

**Note:** The build process automatically applies performance optimizations including:
- Deferred loading of non-critical JavaScript
- Lazy loading of images  
- Pre-compression (gzip) of static assets
- Cache control headers configuration

See [PERFORMANCE_OPTIMIZATIONS.md](PERFORMANCE_OPTIMIZATIONS.md) for details.

## Check the generated page
To check the generated html run:
```
# For Linux, MacOS
firefox _build/index.html
```


## Upload the generated files to the server
To upload the webpage you built to syscop.de, run:
```
# For Linux, MacOS -- your ssh key must be authorized on the server.
make upload
```
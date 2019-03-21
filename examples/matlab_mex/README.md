This folder contains some examples to use the mex-based acados_matlab interface from Matlab or Octave.
This interface uses the shared libraries created using the make command from the main acados folder
```
make acados_c_shared
```

To run the examples, navigate into the selected folder, and there run the command
```
export ACADOS_INSTALL_DIR="<your acados repo dir>"
# from example directory
source env.sh
```
to set the necessary environment variables.
If `ACADOS_INSTALL_DIR` is not speficied, it will be assumed that the examples are runned from the sub-folders of the current folder (i.e. acados main folder is 2 folders up from the current folder).

Afterwards, launch Matlab or Octave from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensures the correct setting of the environment variables.

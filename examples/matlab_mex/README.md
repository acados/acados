This folder contains some examples to use the mex-based acados_matlab interface.
This interface uses the shared libraries created using the make command from the main acados folder
```
make acados_c_shared
```

To run the examples, navigate into the selected folder, and there run the run the command
```
source ./setenv.sh
```
to set the necessary environment variables.
Afterwards, launch matlab from the same shell.

If you want to run the examples in a different folder, please close the current shell and open a new one to repeat the procedure: this ensure the correct setting of the environment variables.

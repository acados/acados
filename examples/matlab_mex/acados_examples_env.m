
example_dir = fileparts(which('acados_examples_env'));

acados_dir = fullfile(example_dir, '..', '..');

matlab_interface_dir = fullfile(acados_dir, 'interfaces', 'acados_matlab');
addpath(matlab_interface_dir);

casadi_dir = fullfile(acados_dir, 'external', 'casadi-matlab');
addpath(casadi_dir);

setenv('ACADOS_INSTALL_DIR', acados_dir);
setenv('ENV_RUN', 'true');
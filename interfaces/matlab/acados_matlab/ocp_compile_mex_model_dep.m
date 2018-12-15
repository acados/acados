function compile_mex_model_dep(opts_struct)

% get acados folder (if set)
acados_folder = getenv('ACADOS_FOLDER');
% default folder
if length(acados_folder) == 0
	acados_folder = '../../../';
end
% set paths
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' acados_folder, 'interfaces'];
acados_lib_path = ['-L' acados_folder, 'lib'];
acados_matlab_lib_path = ['-L' acados_folder, 'interfaces/matlab/acados_matlab/'];

mex_flags = 'GCC=/usr/bin/gcc-4.9';

% get pointers for external functions in model
if (strcmp(opts_struct.sim_solver, 'erk'))
	mex(mex_flags, acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lmodel', 'ocp_expl_ext_fun_create.c');
elseif (strcmp(opts_struct.sim_solver, 'irk'))
	mex(mex_flags, acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lmodel', 'ocp_impl_ext_fun_create.c');
else
	fprintf('\ncodegen_model: sim_solver not supported: %s\n', opts_struct.sim_solver);
	return;
end


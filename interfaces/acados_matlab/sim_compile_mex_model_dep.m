function sim_compile_mex_model_dep(model_struct, opts_struct)

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = [acados_folder, '/interfaces/acados_matlab/'];
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' acados_folder, '/interfaces'];
acados_lib_path = ['-L' acados_folder, '/lib'];
acados_matlab_lib_path = ['-L' acados_folder, '/interfaces/acados_matlab/'];
model_lib_path = ['-L', pwd];

%% select files to compile
mex_files = {};
% dynamics
if (strcmp(opts_struct.method, 'erk'))
	mex_files = {mex_files{:},
		[acados_mex_folder, 'sim_set_ext_fun_expl.c']
		};
elseif (strcmp(opts_struct.method, 'irk'))
	mex_files = {mex_files{:},
		[acados_mex_folder, 'sim_set_ext_fun_impl.c']
		};
else
	fprintf('\ncodegen_model: method not supported: %s\n', opts_struct.method);
end

%% get pointers for external functions in model
for ii=1:length(mex_files)
	mex(mex_flags, 'CFLAGS=\$CFLAGS -std=c99 -fopenmp', acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, model_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lsim_model', mex_files{ii});
end

function sim_compile_mex_model_dep(model_struct, opts_struct)

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

%% select files to compile
mex_files = {};
% dynamics
if (strcmp(opts_struct.method, 'erk'))
	mex_files = {mex_files{:},
		'sim_expl_ext_fun_create.c' % TODO create sim_set_ext_fun_expl
		};
elseif (strcmp(opts_struct.method, 'irk'))
	mex_files = {mex_files{:},
		'sim_impl_ext_fun_create.c' % TODO create sim_set_ext_fun_impl
		};
else
	fprintf('\ncodegen_model: method not supported: %s\n', opts_struct.method);
end

% to avoid warining on R2017a
mex_flags = 'GCC=/usr/bin/gcc-4.9';

%% get pointers for external functions in model
for ii=1:length(mex_files)
	mex(mex_flags, acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lsim_model', mex_files{ii});
end


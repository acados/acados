function ocp_compile_mex_model_dep(model_struct, opts_struct)

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
if (strcmp(opts_struct.sim_method, 'erk'))
	mex_files = {mex_files{:}, ...
		'ocp_set_ext_fun_expl.c'
		};
elseif (strcmp(opts_struct.sim_method, 'irk'))
	mex_files = {mex_files{:}, ...
		'ocp_set_ext_fun_impl.c'
		};
else
	fprintf('\ncodegen_model: sim_method not supported: %s\n', opts_struct.sim_method);
end
% nonlinear constraints
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'nh') && model_struct.nh>0))
	mex_files = {mex_files{:}, ...
		'ocp_set_ext_fun_h.c'
		};
end
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'nh_e') && model_struct.nh_e>0))
	mex_files = {mex_files{:}, ...
		'ocp_set_ext_fun_h_e.c'
		};
end
% nonlinear least squares
if (strcmp(model_struct.cost_type, 'nls'))
	mex_files = {mex_files{:}, ...
		'ocp_set_ext_fun_y.c'
		};
end

% to avoid warining on R2017a
mex_flags = 'GCC=/usr/bin/gcc-4.9';

%% get pointers for external functions in model
for ii=1:length(mex_files)
	mex(mex_flags, acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-locp_model', mex_files{ii});
end


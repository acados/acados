function ocp_compile_mex_model_dep(model_struct, opts_struct)

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = [acados_folder, '/interfaces/acados_matlab/'];
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' acados_folder, '/interfaces'];
external_include = ['-I' acados_folder, '/external'];
blasfeo_include = ['-I' acados_folder, '/external/blasfeo/include'];
acados_lib_path = ['-L' acados_folder, '/lib'];
acados_matlab_lib_path = ['-L' acados_folder, '/interfaces/acados_matlab/'];
model_lib_path = ['-L', pwd];

%% select files to compile
mex_files = {};
% dynamics
if (strcmp(opts_struct.sim_method, 'erk'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_expl.c']
		};
elseif (strcmp(opts_struct.sim_method, 'irk'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_impl.c']
		};
else
	fprintf('\ncodegen_model: sim_method not supported: %s\n', opts_struct.sim_method);
end
% nonlinear constraints
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'nh') && model_struct.nh>0))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_h.c']
		};
end
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'nh_e') && model_struct.nh_e>0))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_h_e.c']
		};
end
% nonlinear least squares
if (strcmp(model_struct.cost_type, 'nonlinear_ls'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_y.c']
		};
end
if (strcmp(model_struct.cost_e_type, 'nonlinear_ls'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_y_e.c']
		};
end
% external cost
if (strcmp(model_struct.cost_type, 'ext_cost'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_ext_cost.c']
		};
end
if (strcmp(model_struct.cost_e_type, 'ext_cost'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_ext_cost_e.c']
		};
end

%% get pointers for external functions in model
for ii=1:length(mex_files)
	disp(['compiling ', mex_files{ii}])
	mex(mex_flags, 'CFLAGS=\$CFLAGS -std=c99 -fopenmp', acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, acados_matlab_lib_path, model_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-locp_model', mex_files{ii});
end


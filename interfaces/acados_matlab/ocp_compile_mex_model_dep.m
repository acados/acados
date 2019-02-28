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
		[acados_mex_folder, 'ocp_set_ext_fun_dyn_expl.c']
		};
elseif (strcmp(opts_struct.sim_method, 'irk'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_dyn_impl.c']
		};
else
	fprintf('\nocp_compile_mex_mode_dep: sim_method not supported: %s\n', opts_struct.sim_method);
end
% nonlinear constraints
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'dim_nh') && model_struct.dim_nh>0))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_constr_h.c']
		};
end
if (strcmp(model_struct.constr_type, 'bgh') && (isfield(model_struct, 'dim_nh_e') && model_struct.dim_nh_e>0))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_constr_h_e.c']
		};
end
% nonlinear least squares
if (strcmp(model_struct.cost_type, 'nonlinear_ls'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_cost_y.c']
		};
end
if (strcmp(model_struct.cost_type_e, 'nonlinear_ls'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_cost_y_e.c']
		};
end
% external cost
if (strcmp(model_struct.cost_type, 'ext_cost'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_cost_ext_cost.c']
		};
end
if (strcmp(model_struct.cost_type_e, 'ext_cost'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'ocp_set_ext_fun_cost_ext_cost_e.c']
		};
end

if is_octave()
	if exist('cflags_octave.txt')==0
		diary 'cflags_octave.txt'
		diary on
		mkoctfile -p CFLAGS
		diary off
		input_file = fopen('cflags_octave.txt', 'r');
		cflags_tmp = fscanf(input_file, '%[^\n]s');
		fclose(input_file);
		cflags_tmp = [cflags_tmp, ' -std=c99 -fopenmp'];
		input_file = fopen('cflags_octave.txt', 'w');
		fprintf(input_file, '%s', cflags_tmp);
		fclose(input_file);
	end
	input_file = fopen('cflags_octave.txt', 'r');
	cflags_tmp = fscanf(input_file, '%[^\n]s');
	fclose(input_file);
	setenv('CFLAGS', cflags_tmp);
end

%% get pointers for external functions in model
for ii=1:length(mex_files)
	disp(['compiling ', mex_files{ii}])
	if is_octave()
%		mkoctfile -p CFLAGS
		mex(acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, acados_matlab_lib_path, model_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-locp_model', mex_files{ii});
	else
		mex(mex_flags, 'CFLAGS=\$CFLAGS -std=c99 -fopenmp', acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, acados_matlab_lib_path, model_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-locp_model', mex_files{ii});
	end
end


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
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'sim_set_ext_fun_dyn_expl.c'] ...
		};
elseif (strcmp(opts_struct.method, 'irk'))
	mex_files = {mex_files{:}, ...
		[acados_mex_folder, 'sim_set_ext_fun_dyn_impl.c'] ...
		};
else
	fprintf('\ncodegen_model: method not supported: %s\n', opts_struct.method);
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
		mex(acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, model_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lsim_model', mex_files{ii});
	else
		mex(mex_flags, 'CFLAGS=\$CFLAGS -std=c99 -fopenmp', acados_include, acados_interfaces_include, acados_lib_path, acados_matlab_lib_path, model_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', '-lsim_model', mex_files{ii});
	end
end

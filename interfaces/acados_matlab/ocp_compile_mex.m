function ocp_compile_mex()

% get acados folder
acados_folder = getenv('ACADOS_INSTALL_DIR');
mex_flags = getenv('ACADOS_MEX_FLAGS');

% set paths
acados_mex_folder = [acados_folder, '/interfaces/acados_matlab/'];
acados_include = ['-I', acados_folder];
acados_interfaces_include = ['-I', acados_folder, '/interfaces'];
external_include = ['-I', acados_folder, '/external'];
blasfeo_include = ['-I', acados_folder, '/external/blasfeo/include'];
acados_lib_path = ['-L', acados_folder, '/lib'];

% compile mex
mex_files = { ...
	[acados_mex_folder, 'ocp_create.c'], ...
	[acados_mex_folder, 'ocp_destroy.c'], ...
	[acados_mex_folder, 'ocp_create_ext_fun.c'], ...
	[acados_mex_folder, 'ocp_destroy_ext_fun.c'], ...
	[acados_mex_folder, 'ocp_solve.c'], ...
	[acados_mex_folder, 'ocp_set.c'], ...
	[acados_mex_folder, 'ocp_get.c'], ...
	[acados_mex_folder, 'ocp_set_model.c'] ...
	} ;


for ii=1:length(mex_files)
	disp(['compiling ', mex_files{ii}])
	mex(mex_flags, 'CFLAGS=\$CFLAGS -std=c99 -fopenmp', acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', mex_files{ii})
end

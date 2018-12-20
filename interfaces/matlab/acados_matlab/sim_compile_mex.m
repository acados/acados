function sim_compile_mex()

% get acados folder
acados_folder = getenv('ACADOS_FOLDER');

% set paths
acados_mex_folder = [acados_folder, '/interfaces/matlab/acados_matlab/'];
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' acados_folder, '/interfaces'];
acados_lib_path = ['-L' acados_folder, '/lib'];

% compile mex
mex_files ={ ...
	[acados_mex_folder, 'sim_create.c'], ...
	[acados_mex_folder, 'sim_destroy.c'], ...
	[acados_mex_folder, 'sim_destroy_ext_fun.c'], ...
	[acados_mex_folder, 'sim_solve.c'], ...
	[acados_mex_folder, 'sim_set.c'], ...
	[acados_mex_folder, 'sim_get.c'], ...
	[acados_mex_folder, 'sim_set_model.c'] ...
	} ;

% to avoid warnings on R2017a
mex_flags = 'GCC=/usr/bin/gcc-4.9';

for ii=1:length(mex_files)
	mex(mex_flags, acados_include, acados_interfaces_include, acados_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', mex_files{ii})
end



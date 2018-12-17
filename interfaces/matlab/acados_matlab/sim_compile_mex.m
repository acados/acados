function sim_compile_mex()

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

% compile mex
mex_files ={'sim_create.c',
	'sim_destroy.c',
	'sim_ext_fun_destroy.c',
	'sim_solve.c',
	'sim_set.c',
	'sim_get.c',
	'sim_set_model.c'} ;

% to avoid warnings on R2017a
mex_flags = 'GCC=/usr/bin/gcc-4.9';

for ii=1:length(mex_files)
	mex(mex_flags, acados_include, acados_interfaces_include, acados_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', mex_files{ii})
end



function compile_mex()

% get acados folder (if set)
acados_folder = getenv('ACADOS_FOLDER');
% default folder
if length(acados_folder) == 0
	acados_folder = '../../../';
end
% set paths
acados_include = ['-I' acados_folder];
acados_interfaces_include = ['-I' acados_folder, 'interfaces'];
external_include = ['-I' acados_folder, 'external'];
blasfeo_include = ['-I' acados_folder, 'external/blasfeo/include'];
acados_lib_path = ['-L' acados_folder, 'lib'];

% compile mex
mex_files ={
	'ocp_create.c',
	'ocp_destroy.c',
	'ocp_ext_fun_destroy.c',
	'ocp_solve.c',
%	'sim_set.c',
	'ocp_get.c',
	'ocp_set_model.c'
	} ;

mex_flags = 'GCC=/usr/bin/gcc-4.9';


for ii=1:length(mex_files)
	mex(mex_flags, acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', mex_files{ii})
end



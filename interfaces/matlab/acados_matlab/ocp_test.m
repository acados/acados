%% test of native matlab interface

%% compile mex files
% mex -v GCC='/usr/bin/gcc-4.9' ... (-v for verbose, GCC=... to change compiler)

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
mex_files ={
%	'sim_create.c',
%	'sim_destroy.c',
%	'sim_ext_fun_destroy.c',
%	'sim_solve.c',
%	'sim_set.c',
%	'sim_get.c',
%	'sim_set_model.c'
	} ;

for ii=1:length(mex_files)
	mex(acados_include, acados_interfaces_include, acados_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', mex_files{ii})
end



%% arguments
codgen_model = 'true';
sim_scheme = 'erk';



%% model
%model_funs = linear_model;
dyn_model = linear_mass_spring_model;
%model_funs = crane_model;

Ts = 0.5;
N = 10;
nx = dyn_model.nx;
nu = dyn_model.nu;
nbx = nx/2;
nbu = nu;



%% acados ocp model
ocp_model = acados_ocp_model();
ocp_model.set('N', N);
ocp_model.set('Ts', Ts);
ocp_model.set('nx', nx);
ocp_model.set('nu', nu);
ocp_model.set('nbx', nbx);
ocp_model.set('nbu', nbu);
ocp_model.model_struct




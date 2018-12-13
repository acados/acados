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
external_include = ['-I' acados_folder, 'external'];
blasfeo_include = ['-I' acados_folder, 'external/blasfeo/include'];
acados_lib_path = ['-L' acados_folder, 'lib'];

% compile mex
mex_files ={
	'ocp_create.c',
	'ocp_destroy.c',
%	'sim_ext_fun_destroy.c',
%	'sim_solve.c',
%	'sim_set.c',
%	'sim_get.c',
%	'sim_set_model.c'
	} ;

for ii=1:length(mex_files)
	mex(acados_include, acados_interfaces_include, external_include, blasfeo_include, acados_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', mex_files{ii})
end



%% arguments
codgen_model = 'true';
nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
qp_solver_N_pcond = 5;
sim_solver = 'erk';
sim_solver_num_stages = 4; % TODO
sim_solver_num_steps = 3; % TODO



%% model
%model_funs = linear_model;
dyn_model = linear_mass_spring_model;
%model_funs = crane_model;

Ts = 0.5;
N = 10;
nx = dyn_model.nx;
nu = dyn_model.nu;
nyl = nu+nx; % number of outputs in lagrange term
nym = nx; % number of outputs in mayer term
nbx = nx/2;
nbu = nu;

Wl = eye(nyl); for ii=1:nu Wl(ii,ii)=2.0; end
Wm = eye(nym);



%% acados ocp model
ocp_model = acados_ocp_model();
ocp_model.set('N', N);
ocp_model.set('Ts', Ts);
ocp_model.set('nx', nx);
ocp_model.set('nu', nu);
ocp_model.set('nyl', nyl);
ocp_model.set('nym', nym);
ocp_model.set('nbx', nbx);
ocp_model.set('nbu', nbu);
ocp_model.set('Wl', Wl);
ocp_model.set('Wm', Wm);
ocp_model.model_struct



% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_N_pcond', qp_solver_N_pcond);
ocp_opts.set('sim_solver', sim_solver);
ocp_opts.set('sim_solver_num_stages', sim_solver_num_stages);
ocp_opts.set('sim_solver_num_steps', sim_solver_num_steps);
ocp_opts.opts_struct



ocp = acados_ocp(ocp_model, ocp_opts);
ocp
ocp.C_ocp



fprintf('\nsuccess!\n\n');


return;

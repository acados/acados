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
param_scheme = 'multiple_shooting_unif_grid';
param_scheme_N = 10;
nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
qp_solver_N_pcond = 5;
sim_solver = 'erk';
sim_solver_num_stages = 4;
sim_solver_num_steps = 3;



%% create model entries
%model_funs = linear_model;
dyn_model = linear_mass_spring_model;
%model_funs = crane_model;

% dims
T = 5.0; % horizon length time
nx = dyn_model.nx;
nu = dyn_model.nu;
nyl = nu+nx; % number of outputs in lagrange term
nym = nx; % number of outputs in mayer term
nbx = nx/2;
nbu = nu;
% cost
Vul = zeros(nu, nyl); for ii=1:nu Vul(ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vxl = zeros(nx, nyl); for ii=1:nx Vxl(ii,nu+ii)=1.0; end % state-to-output matrix in lagrange term
Vxm = zeros(nx, nym); for ii=1:nx Vxm(ii,ii)=1.0; end % state-to-output matrix in mayer term
Wl = eye(nyl); for ii=1:nu Wl(ii,ii)=2.0; end % weight matrix in lagrange term
Wm = eye(nym); % weight matrix in mayer term
yrl = zeros(nyl, 1); % output reference in lagrange term
yrm = zeros(nym, 1); % output reference in mayer term




%% acados ocp model
ocp_model = acados_ocp_model();
% dims
ocp_model.set('T', T);
ocp_model.set('nx', nx);
ocp_model.set('nu', nu);
ocp_model.set('nyl', nyl);
ocp_model.set('nym', nym);
ocp_model.set('nbx', nbx);
ocp_model.set('nbu', nbu);
% cost
ocp_model.set('Vul', Vul);
ocp_model.set('Vxl', Vxl);
ocp_model.set('Vxm', Vxm);
ocp_model.set('Wl', Wl);
ocp_model.set('Wm', Wm);
ocp_model.set('yrl', yrl);
ocp_model.set('yrm', yrm);
% dynamics

ocp_model.model_struct



% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme', param_scheme);
ocp_opts.set('param_scheme_N', param_scheme_N);
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

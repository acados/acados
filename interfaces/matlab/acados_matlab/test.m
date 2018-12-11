%% test of native matlab interface

%mex -v GCC='/usr/bin/gcc-4.9' ../../../lib/libacore.so ../../../lib/libhpipm.so ../../../lib/libblasfeo.so -lm mex_sim.c
%mex libacore.so libhpipm.so libblasfeo.so sim_create.c

acados_folder = getenv('ACADOS_FOLDER');

if length(acados_folder) == 0
	acados_folder = '../../../';
end

include_acados = ['-I' acados_folder];
include_interfaces = ['-I' acados_folder, 'interfaces'];
acados_lib_path = ['-L' acados_folder, 'lib'];

mex_files ={'sim_create.c',
	'sim_destroy.c',
	'sim_solve.c',
	'sim_set.c',
	'sim_get.c',
	'sim_set_model.c'} ;

for ii=1:length(mex_files)
	mex(include_acados, include_interfaces, acados_lib_path, '-lacados_c', '-lacore', '-lhpipm', '-lblasfeo', mex_files{ii})
end


%% model
model_name = 'model';
%sim_model = crane_model_expl(model_name);
sim_model = linear_model(model_name);

nx = sim_model.nx;
nu = sim_model.nu;


%% acados integrator opts
sim_opts = acados_integrator_opts();
sim_opts.set('codgen_model', 'true');
sim_opts.set('num_stages', 4);
sim_opts.set('num_steps', 3);
sim_opts.set('scheme', 'erk');
sim_opts.set('sens_forw', 'true');


%% acados integrator
% create integrator
sim = acados_integrator(sim_model, sim_opts);
% generate model C functions
sim.codegen_model();
% set input
sim.set('T', 0.5);

x0 = ones(sim_model.nx, 1); %x0(1) = 2.0;
tic;
sim.set('x', x0);
time_set_x = toc

u = ones(nu, 1);
sim.set('u', u);

% solve
tic;
sim.solve();
time_solve = toc


% get TODO with return value !!!!!
% xn
xn = zeros(nx, 1);
sim.get('xn', xn);
xn
% S_forw
S_forw = zeros(nx, nx+nu);
sim.get('S_forw', S_forw);
S_forw
% Sx
Sx = zeros(nx, nx);
sim.get('Sx', Sx);
Sx
% Su
Su = zeros(nx, nu);
sim.get('Su', Su);
Su


fprintf('\nsuccess!\n\n');


return;

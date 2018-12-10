%% test of native matlab interface

%mex -v GCC='/usr/bin/gcc-4.9' ../../../lib/libacore.so ../../../lib/libhpipm.so ../../../lib/libblasfeo.so -lm mex_sim.c
%mex libacore.so libhpipm.so libblasfeo.so sim_create.c
mex -I/home/gianluca/acados/ -I/home/gianluca/acados/interfaces -L/home/gianluca/acados/lib -lacados_c -lacore -lhpipm -lblasfeo sim_create.c
mex -I/home/gianluca/acados/ -I/home/gianluca/acados/interfaces -L/home/gianluca/acados/lib -lacados_c -lacore -lhpipm -lblasfeo sim_destroy.c
%mex -I/home/gianluca/acados/ -I/home/gianluca/acados/interfaces -L/home/gianluca/acados/lib -lacados_c -lacore -lhpipm -lblasfeo sim_ext_fun_create.c



%% model
model_name = 'model';
sim_model = crane_model_expl(model_name);


%% acados integrator opts
sim_opts = acados_integrator_opts();
sim_opts.set('codgen_model', 'true');
sim_opts.set('num_stages', 4);
sim_opts.set('num_steps', 3);
sim_opts.set('scheme', 'erk');
sim_opts.set('sens_forw', 'false');


%% acados integrator
% create integrator
sim = acados_integrator(sim_model, sim_opts);
% generate model C functions
sim.codegen_model();


fprintf('\nsuccess!\n\n');


return;

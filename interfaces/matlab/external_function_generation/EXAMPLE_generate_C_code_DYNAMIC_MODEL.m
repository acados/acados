%% EXAMPLE - GENERATE C CODE of dynamic model for acados
path = pwd;
path_sim_codegen = strcat(path,'/sim');
addpath(path_sim_codegen);

model = TEST_simple_model;
% for ERK integrator
generate_c_code_explicit_ode(model);

% for IRK, LIFTED_IRK integrator
generate_c_code_implicit_ode(model);
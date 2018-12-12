%% EXAMPLE - GENERATE C CODE for ocp_nlp acados
clearvars;

path_sim_codegen = strcat(pwd,'/sim');
path_cost_codegen = strcat(pwd,'/cost');

addpath(path_sim_codegen);
addpath(path_cost_codegen);

model = TEST_simple_model;

%% DYNAMICS
% for ERK integrator
generate_c_code_explicit_ode(model);

% for IRK, LIFTED_IRK integrator
generate_c_code_implicit_ode(model);


%% COST
% for NLS cost
cost.nls_expr = [model.x; model.u];
cost.name = model.name

generate_c_code_nls_cost( model, cost );
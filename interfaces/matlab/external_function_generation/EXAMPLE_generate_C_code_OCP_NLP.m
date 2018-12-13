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
% specify cost for intermediate stages
cost.nls_expr = [model.x; model.u];
cost.name = model.name;
generate_c_code_nls_cost( model, cost);

% specify and generate cost for final stage
opts.is_terminal = 1;
cost.nls_expr = model.x + ones(size(model.x));
generate_c_code_nls_cost( model, cost, opts);

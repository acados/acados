%% EXAMPLE - GENERATE C CODE for ocp_nlp acados
clearvars;

path_sim_codegen = strcat(pwd,'/sim');
path_cost_codegen = strcat(pwd,'/cost');

addpath(path_sim_codegen);
addpath(path_cost_codegen);

model = TEST_simple_model;

nx = length(model.x);
nu = length(model.u);

%% DYNAMICS
% for ERK integrator
generate_c_code_explicit_ode(model);

% for IRK, LIFTED_IRK integrator
generate_c_code_implicit_ode(model);


%% COST
%% cost for intermediate stages
opts.is_terminal = 0;
% NLS - nonlinear least squares
cost.nls_expr = [model.x; model.u];
cost.name = model.name;
generate_c_code_nls_cost( model, cost, opts);


% externally provided cost - general cost expression
cost.general_expr = sum(model.x);
for i = 1:nu
    cost.general_expr = cost.general_expr + model.u(i)^4; % sum all
end
generate_c_code_external_cost( model, cost, opts);

%% cost for terminal stage
% NLS - nonlinear least squares
opts.is_terminal = 1;
cost.nls_expr = model.x + ones(size(model.x));
generate_c_code_nls_cost( model, cost, opts);

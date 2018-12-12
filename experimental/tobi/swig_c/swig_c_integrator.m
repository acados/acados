import acados.*

%% Initialize

NREP = 500;

nx = acados.intp();
nu = acados.intp();
ns = acados.intp();
num_steps = acados.intp();
sens_adj = acados.boolp();
nx.assign(4);
nu.assign(1);
ns.assign(4);
num_steps.assign(5);
sens_adj.assign(true);

%NF = nx + nu; % columns of forward seed

T = 0.05;

% sim plan & config
plan = sim_solver_plan;
plan.sim_solver(ERK);
config = sim_config_create(plan);
% sim dims
dims = sim_dims_create(config);

sim_dims_set(config, dims, 'nx', nx);  % in method 'sim_dims_set', argument 4 of type 'int const *'
sim_dims_set(config, dims, 'nu', nu);

opts = sim_opts_create(config, dims);
acados.sim_opts_set(config, opts, 'ns', ns);
acados.sim_opts_set(config, opts, 'num_steps', num_steps);
acados.sim_opts_set(config, opts, 'sens_adj', sens_adj);

in = sim_in_create(config, dims);
out = sim_out_create(config, dims);


expl_vde_for = external_function_casadi;
%acados.external_function_casadi_set_fun(expl_vde_for, acados.vdeFun)
%external_function_casadi_create(expl_vde_for);
%expl_vde_adj = external_function_casadi;

sim_model_set(config, in, 'expl_vde_for', expl_vde_for);
sim_model_set(config, in, 'expl_vde_adj', expl_vde_adj);

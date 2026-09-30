
%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

addpath('../linear_mass_spring_model/');

N = 20;
shooting_nodes = [ linspace(0,1,N/2) linspace(1.1,5,N/2+1) ];

model_name = 'lin_mass';

nlp_solver = 'SQP';
regularize_method = 'CONVEXIFY';
nlp_solver_max_iter = 100;
nlp_solver_ext_qp_res = 1;
qp_solver = 'PARTIAL_CONDENSING_HPIPM';
qp_solver_cond_N = 5;
sim_method = 'DISCRETE';
sim_method_num_stages = 4 * ones(N,1);
sim_method_num_stages(end) = 2;
sim_method_num_steps = 3;

model = get_linear_mass_spring_model();
model.name = model_name;


T = 10.0;
nx = length(model.x);
nu = length(model.u);
x0 = zeros(nx, 1);
x0(1) = 2.5;
x0(2) = 2.5;
lh = [-0.5 * ones(nu, 1); -4 * ones(nx, 1)];
uh = [0.5 * ones(nu, 1); 4 * ones(nx, 1)];
lh_e = -4 * ones(nx, 1);
uh_e = 4 * ones(nx, 1);

ocp = AcadosOcp();
ocp.model = model;
ocp.solver_options.N_horizon = N;
ocp.solver_options.shooting_nodes = shooting_nodes;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = nlp_solver;
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.regularize_method = regularize_method;
ocp.solver_options.nlp_solver_ext_qp_res = nlp_solver_ext_qp_res;
ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter;
ocp.solver_options.qp_solver = qp_solver;
ocp.solver_options.integrator_type = sim_method;
ocp.solver_options.sim_method_num_stages = sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = sim_method_num_steps;
if strcmp(qp_solver, 'PARTIAL_CONDENSING_HPIPM')
	ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
end
ocp.code_gen_options.ext_fun_compile_flags = '';

ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.cost.cost_type_0 = 'EXTERNAL';
ocp.model.cost_expr_ext_cost_0 = 0.5 * model.u' * (2 * eye(nu)) * model.u;
ocp.model.cost_expr_ext_cost = 0.5 * ([model.u; model.x] .* [2*ones(nu,1); ones(nx,1)])' * ...
    ([model.u; model.x] .* [2*ones(nu,1); ones(nx,1)]);
ocp.model.cost_expr_ext_cost_e = 0.5 * model.x' * model.x;

ocp.model.con_h_expr = [model.u; model.x];
ocp.model.con_h_expr_e = model.x;
ocp.constraints.lh = lh;
ocp.constraints.uh = uh;
ocp.constraints.lh_e = lh_e;
ocp.constraints.uh_e = uh_e;
ocp.constraints.x0 = x0;

ocp_solver = AcadosOcpSolver(ocp);


x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);


tic;
ocp_solver.solve();
time_ext = toc;

filename = 'iterate.json';
ocp_solver.store_iterate(filename, true);
ocp_solver.load_iterate(filename);
delete(filename)

filename = 'qp.json';
ocp_solver.dump_last_qp_to_json(filename, false, 'Matlab')
delete(filename)

filename = 'qp_c_backend.json';
ocp_solver.dump_last_qp_to_json(filename, false)
delete(filename)

qp_diagnostics_result = ocp_solver.qp_diagnostics();

u = ocp_solver.get('u');
x = ocp_solver.get('x');

status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

ocp_solver.print('stat');

if status~=0
    error('ocp_nlp solver returned status nonzero');
elseif sqp_iter > 2
    error('ocp can be solved in 2 iterations!');
else
	fprintf(['\ntest_ocp_linear_mass_spring: success with sim method ', ...
        sim_method, ' !\n']);
end

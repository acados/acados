%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

import casadi.*
addpath('../pendulum_on_cart_model/');


N = 20;
T = 1;
x0 = [0; pi; 0; 0];

nlp_solver = 'SQP';
qp_solver = 'PARTIAL_CONDENSING_OSQP';
qp_solver_cond_N = 5;
sim_method = 'ERK';

model = get_pendulum_on_cart_model();
nx = length(model.x);
nu = length(model.u);
model.name = 'pendulum';
model.cost_expr_ext_cost = 0.5 * (model.x' * diag([1e3, 1e3, 1e-2, 1e-2]) * model.x ...
    + model.u' * 1e-2 * model.u);
model.cost_expr_ext_cost_e = 0.5 * model.x' * diag([1e3, 1e3, 1e-2, 1e-2]) * model.x;
model.con_h_expr = model.u;
model.con_h_expr_0 = model.u;

ocp = AcadosOcp();
ocp.model = model;
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = nlp_solver;
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.integrator_type = sim_method;
ocp.solver_options.qp_solver = qp_solver;
ocp.solver_options.qp_solver_iter_max = 2000;
ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
ocp.code_gen_options.ext_fun_compile_flags = '';
ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_0 = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.model.cost_expr_ext_cost_0 = 0.5 * model.u' * 1e-2 * model.u;

U_max = 80;
ocp.constraints.lh = -U_max;
ocp.constraints.uh = U_max;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh_0 = U_max;
ocp.constraints.x0 = x0;

ocp_solver = AcadosOcpSolver(ocp);

x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);

ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

ocp_solver.solve();
utraj = ocp_solver.get('u');
xtraj = ocp_solver.get('x');

status = ocp_solver.get('status');
ocp_solver.print('stat');

if status == 0
    disp('test_ocp_OSQP: success!');
else
    error(['test_ocp_OSQP: Failed with status ', num2str(status)]);
end

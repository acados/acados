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
qp_solver = 'PARTIAL_CONDENSING_QPDUNES';
qp_solver_cond_N = 5;
sim_method = 'ERK';

model = get_pendulum_on_cart_model();
nx = length(model.x);
nu = length(model.u);
model.name = 'pendulum';
model.con_h_expr = model.u;
model.con_h_expr_0 = model.u;

cost_type = 'LINEAR_LS';
ny = nx + nu;
ny_e = nx;

Vu = zeros(ny, nu);
for ii = 1:nu
    Vu(ii,ii) = 1.0;
end
Vx = zeros(ny, nx);
for ii = 1:nx
    Vx(nu+ii,ii) = 1.0;
end
Vx_e = eye(ny_e, nx);
W = diag([1, 100, 100, 1, 1]);
W_e = W(nu+1:nu+nx, nu+1:nu+nx);
yr = zeros(ny, 1);
yr_e = zeros(ny_e, 1);

ocp = AcadosOcp();
ocp.model = model;
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = nlp_solver;
ocp.solver_options.integrator_type = sim_method;
ocp.solver_options.qp_solver = qp_solver;
ocp.solver_options.qp_solver_iter_max = 2000;
ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
ocp.code_gen_options.ext_fun_compile_flags = '';
ocp.cost.cost_type = cost_type;
ocp.cost.cost_type_e = cost_type;
ocp.cost.Vu = Vu;
ocp.cost.Vx = Vx;
ocp.cost.Vx_e = Vx_e;
ocp.cost.W = W;
ocp.cost.W_e = W_e;
ocp.cost.yref = yr;
ocp.cost.yref_e = yr_e;

U_max = 80;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh_0 = U_max;
ocp.constraints.lh = -U_max;
ocp.constraints.uh = U_max;
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
    disp('test_ocp_qpDUNES: success!');
else
    error(['test_ocp_qpDUNES: Failed with status ', num2str(status)]);
end

%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

import casadi.*

modelFunction = Function.load('modelFunction');

nStatesAndInputs = 8;
nStates          = 6;
Vx               = [1 0 0 0 0 0 0 0;
                    0 1 0 0 0 0 0 0];
ref              = [0.05;1];
Weights          = diag([400, 130]);
StatesAndInputs0 = [-0.0061;
                     1.0056;
                     0.0409;
                     0.6621;
                     0.0184;
                     0.2067;
                     0.0500;
                     0.0500];

x = SX.sym('Opts',nStatesAndInputs,1);

xdot = modelFunction(x);

J = (Vx*x-ref)'*Weights*(Vx*x-ref);

constr = [xdot;x(7:8)];
LB     = [zeros(nStates,1);0;0];
UB     = [zeros(nStates,1);1;1];

optionsIPOPT = struct('ipopt', struct('max_iter', 500));
prob         = struct('f', J, 'x', x, 'g', constr);
IPOPTFun     = nlpsol('solver', 'ipopt', prob, optionsIPOPT);

sol = IPOPTFun('x0', StatesAndInputs0,'lbg', LB, 'ubg', UB);

x_ipopt = full(sol.x);

N = 1;
T = 1;

nlp_solver = 'SQP';
qp_solver = 'FULL_CONDENSING_HPIPM';
sim_method = 'DISCRETE';

model = AcadosModel();
model.name = 'target_selector';
model.x = x;
model.disc_dyn_expr = x;
model.cost_expr_ext_cost_0 = J;
model.cost_expr_ext_cost = J;
model.cost_expr_ext_cost_e = SX.zeros(1);
model.con_h_expr_e = xdot;

ocp = AcadosOcp();
ocp.model = model;
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = nlp_solver;
ocp.solver_options.integrator_type = sim_method;
ocp.solver_options.qp_solver = qp_solver;
ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_0 = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.constraints.lh_e = zeros(size(xdot));
ocp.constraints.uh_e = zeros(size(xdot));

ocp.simulink_opts = AcadosOcpSimulinkOptions();
ocp.simulink_opts.inputs.x_init = 1;
ocp.simulink_opts.outputs.u0 = 0;
ocp.simulink_opts.outputs.sqp_iter = 0;
ocp.simulink_opts.outputs.CPU_time = 0;
ocp.simulink_opts.outputs.x1 = 0;

ocp_solver = AcadosOcpSolver(ocp);

x_traj_init = repmat(StatesAndInputs0, 1, 2);

ocp_solver.set('init_x', x_traj_init);

ocp_solver.solve();
disp(['acados ocp solver returned status ', ocp_solver.get('status')]);
ocp_solver.print('stat');

x_acados = ocp_solver.get('x', 0);

diff_acados_ipopt = norm(x_acados-x_ipopt)


tol = 1e-6;
if any([diff_acados_ipopt] > tol)
    disp(['diff_acados_ipopt', diff_acados_ipopt'])
    error(['test_target_selector: solution of templated MEX and original MEX and IPOPT',...
         ' differ too much. Should be < tol = ' num2str(tol)]);
end

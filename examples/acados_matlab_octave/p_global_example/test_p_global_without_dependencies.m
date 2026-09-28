%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

import casadi.*
%
check_acados_requirements()

%% solver settings
N = 1; % number of discretization steps
T = 1; % [s] prediction horizon length

%% model dynamics
model = AcadosModel();
model.x = SX.sym('x', 1);
model.p_global = SX.sym('p_global', 500);
model.name = 'test_model';
model.disc_dyn_expr = model.x;

%% OCP formulation object
ocp = AcadosOcp();
ocp.model = model;

% path cost term
ocp.cost.cost_type = 'EXTERNAL';
ocp.model.cost_expr_ext_cost = (model.x-1)^2;

% terminal cost term
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.model.cost_expr_ext_cost_e = 0;

%% define constraints
% only bound on u on initial stage and path
ocp.constraints.lbx_0 = 0;
ocp.constraints.ubx_0 = 0.5;
ocp.constraints.idxbx_0 = 0;

ocp.p_global_values = ones(500, 1);

% define solver options
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.integrator_type = 'DISCRETE';
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';

% create solver
ocp_solver = AcadosOcpSolver(ocp);

%% call ocp solver
% solve
ocp_solver.solve();
% get solution
utraj = ocp_solver.get('u');
xtraj = ocp_solver.get('x');

disp('xtraj')
disp(xtraj)

if abs(xtraj - 0.5) > 1e-6
    error('acados returned wrong solution');
end

status = ocp_solver.get('status'); % 0 - success
ocp_solver.print('stat')

if status ~= 0
    error('acados returned status %d', status);
end

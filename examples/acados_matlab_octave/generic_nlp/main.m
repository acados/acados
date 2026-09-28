%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

clear all; clc; close all
check_acados_requirements()
import casadi.*

%% generic NLP formulation
% min. f(x,p)
% s.t. glb <= g(x,p) <= gub

x = casadi.SX.sym('x',2);
p = casadi.SX.sym('p',2);
f = p(1) * (100*(x(2)-x(1)^2)^2 + (x(1)-1)^2);
g = x(1)^2 + x(2)^2 - 2*p(2);
glb = 0;
gub = 0;

%% acados model object
model = AcadosModel();
model.name = 'generic_nlp';
model.x = x;
model.p = p;

%% acados ocp formulation
ocp = AcadosOcp();
ocp.name = 'nlp_solver';
ocp.model = model;

% (terminal) cost
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.model.cost_expr_ext_cost_e = f;

% (terminal) constraints
ocp.model.con_h_expr_e = g;
ocp.constraints.lh_e = glb;
ocp.constraints.uh_e = gub;

% initial parameter values
ocp.parameter_values = zeros(length(model.p),1);

%% solver options
ocp.solver_options.N_horizon = 0;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';

%% create the solver
ocp_solver = AcadosOcpSolver(ocp);

%% solve the NLP
% initial guess
init_x = [2.5; 3.0];
ocp_solver.set('init_x', init_x);

% set the parameters
p_value = [1;1];
ocp_solver.set('p', p_value, 0, ocp.solver_options.N_horizon+1); % set parameter values for all stages at once

% solve and time
tic
ocp_solver.solve();
time_external = toc;
% internal timing
total_time = ocp_solver.get('time_tot');


% check status
status = ocp_solver.get('status');
if status ~= 0
    warning(['solver failed with status ',num2str(status)]);
end

% display results
x_opt = ocp_solver.get('x', 0);
disp('Optimal solution:')  % should be [1;1] for p = [1;1]
disp(x_opt)
disp(['Total time (internal): ', num2str(1e3*total_time), ' ms'])
disp(['Total time (external): ', num2str(1e3*time_external), ' ms'])

% compare with the expected solution
if all(p_value == [1;1])
    assert(sum(x_opt-[1;1]) < 1e-6, 'solution does not match the expected one')
end

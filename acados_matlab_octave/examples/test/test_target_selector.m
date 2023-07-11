%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

% simple example of a target selector (N=1, goal: find steady state)
% model of a motor air path
% author: Severin Hänggi (& Jonathan Frey)

import casadi.*

% load model
modelFunction = Function.load('modelFunction');

% Define parameters
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

% Define symbolic optimization variables
x = SX.sym('Opts',nStatesAndInputs,1);

% get symbolic state derivatives (equal to equality constraints)
xdot = modelFunction(x);

% Derive symbolic cost
J = (Vx*x-ref)'*Weights*(Vx*x-ref);

% Derive Symbolic Constraints (state derivatives and input)
constr = [xdot;x(7:8)];
LB     = [zeros(nStates,1);0;0];
UB     = [zeros(nStates,1);1;1];

% Set Up IPOPT Prob
optionsIPOPT = struct('ipopt', struct('max_iter', 500));
prob         = struct('f', J, 'x', x, 'g', constr);
IPOPTFun     = nlpsol('solver', 'ipopt', prob, optionsIPOPT);

% Run IPOPT Problem
sol = IPOPTFun('x0', StatesAndInputs0,'lbg', LB, 'ubg', UB);

x_ipopt = full(sol.x);
% sol.x = [0.05, 1, 0.0634343, 0.667782, -0.00725671, 0.269601, 0.334563,
% 0.155177];

%% acados
nx = nStatesAndInputs;

N = 1;
T = 1; % time horizon length - no meaning here, might be relevant for sampling time of Simulink block

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'full_condensing_hpipm';
sim_method = 'discrete'; % erk, irk, irk_gnsf

%% model to create the solver
ocp_model = acados_ocp_model();
model_name = 'target_selector';

%% acados ocp model
ocp_model.set('name', model_name);
ocp_model.set('T', 1);

% symbolics
ocp_model.set('sym_x', x);

% cost
ocp_model.set('cost_type', 'ext_cost'); % ext_cost, auto
ocp_model.set('cost_type_e', 'auto'); % ext_cost, auto
ocp_model.set('cost_expr_ext_cost', J);
ocp_model.set('cost_expr_ext_cost_e', SX.zeros(1));

% % dynamics
ocp_model.set('dyn_type', 'discrete');
ocp_model.set('dyn_expr_phi', x);
% x_{k+1} = x_{k}

% constraints
% impose nonlinear constraint at terminal node.
% Note: currently: nonlinear constr h is per default not enforced at
% initial shooting node, because typically (in optimal control) there are
% initial state constraints
ocp_model.set('constr_type_e', 'auto');
ocp_model.set('constr_expr_h_e', xdot);
ocp_model.set('constr_lh_e', zeros(size(xdot))); % lower bound on h
ocp_model.set('constr_uh_e', zeros(size(xdot)));  % upper bound on h

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('qp_solver', qp_solver);


%% test with simulink options
simulink_opts = get_acados_simulink_opts;
simulink_opts.inputs.x_init = 1;
simulink_opts.outputs.u0 = 0;

simulink_opts.outputs.sqp_iter = 0;
simulink_opts.outputs.CPU_time = 0;
simulink_opts.outputs.x1 = 0;

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts, simulink_opts);

% initialize
% 2 working initializations
eps = 1e-1;
% x_traj_init = repmat(x_ipopt, 1, 2) + eps * ones(nx, N+1);
x_traj_init = repmat(StatesAndInputs0, 1, 2);

%% call ocp solver
% set trajectory initialization
ocp.set('init_x', x_traj_init);
% ocp.set('init_pi', zeros(nx, N))

ocp.solve();
disp(['acados ocp solver returned status ', ocp.get('status')]); % 0 - success
ocp.print('stat')

% get solution
x_acados = ocp.get('x', 0);

% [x_acados, x_ipopt]
diff_acados_ipopt = norm(x_acados-x_ipopt)


tol = 1e-6;
if any([diff_acados_ipopt] > tol)
    disp(['diff_acados_ipopt', diff_acados_ipopt'])
    error(['test_target_selector: solution of templated MEX and original MEX and IPOPT',...
         ' differ too much. Should be < tol = ' num2str(tol)]);
end

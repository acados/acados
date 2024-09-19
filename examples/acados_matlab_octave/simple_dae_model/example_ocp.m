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


check_acados_requirements()

%% solver settings
N = 20; % number of discretization steps
T = 1.; % [s] prediction horizon length
h = T/N;
x0 = [3; -1.8]; % initial state

%% model dynamics
model = simple_dae_model();

nx = length(model.x);
nu = length(model.u);
nz = length(model.z);
ny = nx + nu;

%% OCP formulation object
ocp = AcadosOcp();
ocp.model = model;

%% Options
ocp.solver_options.compile_interface = 'AUTO';

ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;

ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
ocp.solver_options.nlp_solver_max_iter = 10;
ocp.solver_options.nlp_solver_tol_stat = 1e-12;
ocp.solver_options.nlp_solver_tol_eq   = 1e-12;
ocp.solver_options.nlp_solver_tol_ineq = 1e-12;
ocp.solver_options.nlp_solver_tol_comp = 1e-12;
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.nlp_solver_step_length = 0.7;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5;
ocp.solver_options.qp_solver_warm_start = 0;
ocp.solver_options.qp_solver_cond_ric_alg = 1; % 0: dont factorize hessian in the condensing; 1: factorize
ocp.solver_options.qp_solver_ric_alg = 1; % HPIPM specific
ocp.solver_options.integrator_type = 'IRK';
ocp.solver_options.sim_method_num_stages = 6; % scalar or vector of size N_horizon;
ocp.solver_options.sim_method_num_steps = 4; % scalar or vector of size N_horizon;
ocp.solver_options.sim_method_newton_iter = 3; % scalar or vector of size N_horizon;

% selectors for example variants
constr_variant = 1; % 0: x bounds; 1: z bounds
cost_variant = 2; % 0: ls on u,x; 1: ls on u,z; 2: nls on u, z

% Cost
Wu = 1e-3*eye(nu);
Wx = 1e1*eye(nx);
W = [Wu, zeros(nu,nx); zeros(nx,nu), Wx];
Vu = [eye(nu); zeros(nx,nu)];
Vx = [zeros(nu,nx); eye(nx)];
Vx_e = eye(nx);

if cost_variant == 0
    ocp.cost.cost_type = 'LINEAR_LS';
    ocp.cost.Vu = Vu;
    ocp.cost.Vx = Vx;
    ocp.cost.Vz = zeros(ny, nz);
elseif cost_variant == 1
    ocp.cost.cost_type = 'LINEAR_LS';
    ocp.cost.Vu = Vu;
    ocp.cost.Vz = Vx;
    ocp.cost.Vx = zeros(ny, nx);
else
    ocp.cost.cost_type = 'NONLINEAR_LS';
    ocp.model.cost_y_expr = vertcat(model.u, model.z);
end

ocp.cost.cost_type_e = 'LINEAR_LS';
ocp.cost.Vx_e = Vx_e;
ocp.cost.W = W;
ocp.cost.W_e = Wx;

% Constraints
lb = [2; -2];
ub = [4;  2];
ocp.constraints.x0 = x0;

if constr_variant == 0
    ocp.constraints.idxbx = 0:(nx-1);
    ocp.constraints.lbx = lb;
    ocp.constraints.ubx = ub;
else
    ocp.model.con_h_expr = model.z;
    ocp.constraints.lh = lb;
    ocp.constraints.uh = ub;
    ocp.model.con_h_expr_e = model.x;
    ocp.constraints.lh_e = lb;
    ocp.constraints.uh_e = ub;
end

%% Create solver and solve
ocp_solver = AcadosOcpSolver(ocp);

ocp_solver.solve();

stat = ocp_solver.get('stat');

ocp_solver.print('stat');

status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
sqp_time = ocp_solver.get('time_tot');

% get solution for initialization of next NLP
x_traj = ocp_solver.get('x');
u_traj = ocp_solver.get('u');
pi_traj = ocp_solver.get('pi');
z_traj = ocp_solver.get('z');

diff_x_z = x_traj(:,1:N) - z_traj;

max_diff_x_z = max(max(abs(diff_x_z)));
test_tol = 1e-14;
if max_diff_x_z > test_tol
    error(['test_ocp_simple_dae: diff_x_z > ' num2str(test_tol), ' is ' num2str(max_diff_x_z)]);
end

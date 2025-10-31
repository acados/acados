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


clear all; clc;
check_acados_requirements()

%% acados ocp model
model = get_linear_mass_spring_model();
nx = length(model.x);
nu = length(model.u);

%% set up OCP
ocp = AcadosOcp();
ocp.model = model;

T = 10.0; % horizon length time
N = 20;
sim_method = 'DISCRETE';
shooting_nodes = linspace(0,10,N+1);
tol = 1e-8;
casadi_dynamics = 0; % 0=generic, 1=casadi
casadi_cost = 1; % 0=generic, 1=casadi

ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;
ocp.solver_options.nlp_solver_type = 'SQP'; % 'SQP_RTI'
ocp.solver_options.hessian_approx = 'EXACT'; % 'EXACT', 'GAUSS_NEWTON'
ocp.solver_options.regularize_method = 'CONVEXIFY';% NO_REGULARIZE, PROJECT, PROOJECT_REDUC_HESS, MIRROR, CONVEXIFY
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.nlp_solver_max_iter = 100;
ocp.solver_options.nlp_solver_tol_stat = tol;
ocp.solver_options.nlp_solver_tol_eq = tol;
ocp.solver_options.nlp_solver_tol_ineq = tol;
ocp.solver_options.nlp_solver_tol_comp = tol;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5; % for partial condensing
ocp.solver_options.integrator_type = 'DISCRETE';

% cost
ocp.cost.cost_type_0 = 'EXTERNAL';
ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';

sym_x = model.x;
sym_u = model.u;
yr_u = zeros(nu, 1);
yr_x = zeros(nx, 1);
dWu = 2*ones(nu, 1);
dWx = ones(nx, 1);

ymyr_0 = sym_u - yr_u;
ymyr = [sym_u; sym_x] - [yr_u; yr_x];
ymyr_e = sym_x - yr_x;

cost_expr_ext_cost_0 = 0.5 * ymyr_0' * (dWu .* ymyr_0);
cost_expr_ext_cost = 0.5 * ymyr' * ([dWu; dWx] .* ymyr);
cost_expr_ext_cost_e = 0.5 * ymyr_e' * (dWx .* ymyr_e);

if (casadi_cost == 0)
    ocp.cost.cost_ext_fun_type_0 = 'generic';
    ocp.cost.cost_source_ext_cost_0 = 'generic_ext_cost.c';
    ocp.cost.cost_function_ext_cost_0 =  'ext_cost';
    
    ocp.cost.cost_ext_fun_type = 'generic';
    ocp.cost.cost_source_ext_cost = 'generic_ext_cost.c';
    ocp.cost.cost_function_ext_cost =  'ext_cost';
    
    ocp.cost.cost_ext_fun_type_e = 'generic';
    ocp.cost.cost_source_ext_cost_e = 'generic_ext_cost.c';
    ocp.cost.cost_function_ext_cost_e =  'ext_costN';
else
    ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
    ocp.model.cost_expr_ext_cost = cost_expr_ext_cost;
    ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
end

% dynamic
if (casadi_dynamics == 0)
    ocp.model.dyn_ext_fun_type = 'generic';
    ocp.model.dyn_generic_source = 'generic_disc_dyn.c';
    ocp.model.dyn_disc_fun = 'disc_dyn_fun';
    ocp.model.dyn_disc_fun_jac = 'disc_dyn_fun_jac';
    ocp.model.dyn_disc_fun_jac_hess = 'disc_dyn_fun_jac_hess';
else
    ocp.model.dyn_ext_fun_type = 'casadi';
end

% constraints
constr_expr_h_0 = sym_u;
constr_expr_h = [sym_u; sym_x];
constr_expr_h_e = sym_x;
model.con_h_expr_0 = constr_expr_h_0;
model.con_h_expr = constr_expr_h;
model.con_h_expr_e = constr_expr_h_e;

x0 = zeros(nx, 1); x0(1)=2.5; x0(2)=2.5;
ocp.constraints.x0 = x0;

lh_0 = -0.5 * ones(nu, 1);
uh_0 = 0.5 * ones(nu, 1);
ocp.constraints.lh_0 = lh_0;
ocp.constraints.uh_0 = uh_0;

lh = - [ 0.5 * ones(nu, 1); 4.0 * ones(nx, 1)];
uh = + [ 0.5 * ones(nu, 1); 4.0 * ones(nx, 1)];
ocp.constraints.lh = lh;
ocp.constraints.uh = uh;

lh_e = -4.0 * ones(nx, 1);
uh_e = 4.0 * ones(nx, 1);
ocp.constraints.lh_e = lh_e;
ocp.constraints.uh_e = uh_e;
%% acados ocp solver
% create ocp
ocp_solver = AcadosOcpSolver(ocp);

% initial state
ocp_solver.set('constr_x0', x0);

% set trajectory initialization
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% solve
tic;
ocp_solver.solve();
time_ext = toc;

% get solution
utraj = ocp_solver.get('u');
xtraj = ocp_solver.get('x');

% get info
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_reg = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

% print statistics
ocp_solver.print('stat')

if status~=0
    error('ocp_nlp solver returned status nonzero');
elseif sqp_iter > 2
    error('ocp can be solved in 2 iterations!');
else
	fprintf(['\ntest_ocp_linear_mass_spring: success with sim method ', ...
        sim_method, ' !\n']);
end

% % plot result
% figure()
% subplot(2, 1, 1)
% plot(0:N, xtraj);
% title('trajectory')
% ylabel('x')
% subplot(2, 1, 2)
% plot(1:N, utraj);
% ylabel('u')
% xlabel('sample')

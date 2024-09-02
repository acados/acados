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



% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)


clear all; clc;
check_acados_requirements()
addpath('../linear_mass_spring_model/');

%% arguments
N = 20;
tol = 1e-8;
shooting_nodes = linspace(0,10,N+1);

model_name = 'lin_mass';

nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
% nlp_solver_exact_hessian = 'false';
nlp_solver_exact_hessian = 'true';
% regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'mirror';
regularize_method = 'convexify';
nlp_solver_max_iter = 100;
nlp_solver_ext_qp_res = 1;
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
qp_solver_cond_N = 5;

sim_method = 'discrete';
dyn_type = 'discrete';
cost_type = 'ext_cost';

%% create model entries
model = linear_mass_spring_model();

% dims
T = 10.0; % horizon length time
nx = model.nx;
nu = model.nu;

% constraints
x0 = zeros(nx, 1); x0(1)=2.5; x0(2)=2.5;

lh_0 = -0.5 * ones(nu, 1);
uh_0 = 0.5 * ones(nu, 1);
lh = - [ 0.5 * ones(nu, 1); 4.0 * ones(nx, 1)];
uh = + [ 0.5 * ones(nu, 1); 4.0 * ones(nx, 1)];
lh_e = -4.0 * ones(nx, 1);
uh_e = 4.0 * ones(nx, 1);


%% acados ocp model
casadi_dynamics = 0; % 0=generic, 1=casadi
casadi_cost = 1; % 0=generic, 1=casadi

ocp_model = acados_ocp_model();
ocp_model.set('name', model_name);
ocp_model.set('T', T);

% symbolics
ocp_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
	ocp_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_xdot')
	ocp_model.set('sym_xdot', model.sym_xdot);
end

% cost
ocp_model.set('cost_type_0', cost_type);
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);

% dynamics
ocp_model.set('dyn_type', 'discrete');

if (casadi_dynamics == 0)
    % Generic dynamics
    ocp_model.set('dyn_ext_fun_type', 'generic');
    ocp_model.set('dyn_generic_source', 'generic_disc_dyn.c');
    ocp_model.set('dyn_disc_fun', 'disc_dyn_fun');
    ocp_model.set('dyn_disc_fun_jac', 'disc_dyn_fun_jac');
    ocp_model.set('dyn_disc_fun_jac_hess', 'disc_dyn_fun_jac_hess'); % only needed for exact hessian
else
    % dynamics expression
    ocp_model.set('dyn_expr_phi', model.dyn_expr_phi);
end

if (casadi_cost == 0)
    % Generic initial cost
    ocp_model.set('cost_ext_fun_type_0', 'generic');
    ocp_model.set('cost_source_ext_cost_0', 'generic_ext_cost.c');
    ocp_model.set('cost_function_ext_cost_0', 'ext_cost');
    % Generic stage cost
    ocp_model.set('cost_ext_fun_type', 'generic');
    ocp_model.set('cost_source_ext_cost', 'generic_ext_cost.c');
    ocp_model.set('cost_function_ext_cost', 'ext_cost');
    % Generic terminal cost
    ocp_model.set('cost_ext_fun_type_e', 'generic');
    ocp_model.set('cost_source_ext_cost_e', 'generic_ext_cost.c');
    ocp_model.set('cost_function_ext_cost_e', 'ext_costN');
else
    % cost expression
    ocp_model.set('cost_expr_ext_cost_0', model.cost_expr_ext_cost_0);
    ocp_model.set('cost_expr_ext_cost', model.cost_expr_ext_cost);
    ocp_model.set('cost_expr_ext_cost_e', model.cost_expr_ext_cost_e);
end

% constraints
ocp_model.set('constr_x0', x0);
ocp_model.set('constr_expr_h_0', model.constr_expr_h_0);
ocp_model.set('constr_lh_0', lh_0);
ocp_model.set('constr_uh_0', uh_0);
ocp_model.set('constr_expr_h', model.constr_expr_h);
ocp_model.set('constr_lh', lh);
ocp_model.set('constr_uh', uh);
ocp_model.set('constr_expr_h_e', model.constr_expr_h_e);
ocp_model.set('constr_lh_e', lh_e);
ocp_model.set('constr_uh_e', uh_e);


%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
if (exist('shooting_nodes', 'var'))
	ocp_opts.set('shooting_nodes', shooting_nodes);
end
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);
if (strcmp(nlp_solver, 'sqp')) % not available for sqp_rti
    ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
    ocp_opts.set('nlp_solver_tol_stat', tol);
    ocp_opts.set('nlp_solver_tol_eq', tol);
    ocp_opts.set('nlp_solver_tol_ineq', tol);
    ocp_opts.set('nlp_solver_tol_comp', tol);
end
ocp_opts.set('qp_solver', qp_solver);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
	ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
end

ocp_opts.set('sim_method', sim_method);

%% acados ocp
% create ocp
ocp_solver = acados_ocp(ocp_model, ocp_opts);

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

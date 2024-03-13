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

%% test of native matlab interface
clear all

addpath('../linear_mass_spring_model/');

%% arguments
N = 20;
tol = 1e-8;
shooting_nodes = [ linspace(0,10,N+1) ];

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
    ocp_model.set('dyn_expr_phi', model.expr_phi);
end

if (casadi_cost == 0)
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
    ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
    ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
end

% constraints
ocp_model.set('constr_x0', x0);
ocp_model.set('constr_expr_h', model.expr_h);
ocp_model.set('constr_lh', lh);
ocp_model.set('constr_uh', uh);
ocp_model.set('constr_expr_h_e', model.expr_h_e);
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
ocp = acados_ocp(ocp_model, ocp_opts);

% initial state
ocp.set('constr_x0', x0);

% set trajectory initialization
x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);
ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);


% solve
tic;
ocp.solve();
time_ext = toc;

% get solution
utraj = ocp.get('u');
xtraj = ocp.get('x');

% get info
status = ocp.get('status');
sqp_iter = ocp.get('sqp_iter');
time_tot = ocp.get('time_tot');
time_lin = ocp.get('time_lin');
time_reg = ocp.get('time_reg');
time_qp_sol = ocp.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

% print statistics
ocp.print('stat')

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

%% test templated ocp solver
disp('testing templated solver');
ocp.generate_c_code;
cd c_generated_code/

t_ocp = lin_mass_mex_solver;

% initial state
t_ocp.set('constr_x0', x0);
% t_ocp.set('print_level', print_level)

% set trajectory initialization
t_ocp.set('init_x', x_traj_init);
t_ocp.set('init_u', u_traj_init);

t_ocp.solve();
xt_traj = t_ocp.get('x');
ut_traj = t_ocp.get('u');
status = t_ocp.get('status');

if status~=0
    error('test_template_pendulum_ocp: solution of templated MEX failed!');
else
    fprintf('\ntest_template_pendulum_ocp: templated MEX success!\n');
end

error_X_mex_vs_mex_template = max(max(abs(xt_traj - xtraj)))
error_U_mex_vs_mex_template = max(max(abs(ut_traj - utraj)))

t_ocp.print('stat')


tol_check = 1e-6;

if any([error_X_mex_vs_mex_template, error_U_mex_vs_mex_template] > tol_check)
    error(['test_template_pendulum_exact_hess: solution of templated MEX and original MEX',...
         ' differ too much. Should be < tol = ' num2str(tol_check)]);
end

cost_native_mex = ocp.get_cost()
cost_template_mex = t_ocp.get_cost()

if any(abs(cost_native_mex - cost_template_mex) > tol_check)
    error(['test_template_pendulum_exact_hess: cost value of templated MEX and original MEX',...
         ' differ too much. Should be < tol = ' num2str(tol_check)]);
end

clear all
cd ..

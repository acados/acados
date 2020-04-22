%
% Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
% Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
% Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
% Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

addpath('../pendulum_on_cart_model')

check_acados_requirements()

%% discretization
N = 100;
h = 0.01;
T = N*h; % time horizon length

nlp_solver = 'sqp'; % sqp, sqp_rti
nlp_solver_exact_hessian = 'true';
regularize_method = 'project_reduc_hess';
     % no_regularize, project, project_reduc_hess, mirror, convexify
nlp_solver_max_iter = 100;
tol = 1e-8;
qp_solver = 'partial_condensing_hpipm';
    % full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases
qp_solver_cond_N = 5; % for partial condensing
qp_solver_cond_ric_alg = 0;
qp_solver_ric_alg = 0;
qp_solver_warm_start = 1; % 0: cold, 1: warm, 2: hot
qp_solver_iter_max = 100;
sim_method_num_stages = 4;
sim_method_num_steps = 1;


%% model dynamics
model = pendulum_on_cart_model;

%% model to create the solver
ocp_model = acados_ocp_model();

%% dimensions
nx = model.nx;
nu = model.nu;

nbx = 0;
nbu = 0;
nh = nu;
nh_e = 0;

%% cost formulation
cost_formulation = 3;
switch cost_formulation
    case 1
        cost_type = 'linear_ls';
    case 2
        cost_type = 'ext_cost';
    otherwise
        cost_type = 'auto';
end

%% integrator type
integrator = 2;
switch integrator
    case 1
        sim_method = 'erk';
    case 2
        sim_method = 'irk';
    otherwise
        sim_method = 'irk_gnsf';
end
model_name = ['pendulum_' num2str(cost_formulation) num2str(integrator)];

%% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);
if strcmp( cost_type, 'linear_ls' )
    ny = nu+nx; % number of outputs in lagrange term
    ny_e = nx; % number of outputs in mayer term
    Vu = zeros(ny, nu); for ii=1:nu Vu(ii,ii)=1.0; end % input-to-output matrix in lagrange term
    Vx = zeros(ny, nx); for ii=1:nx Vx(nu+ii,ii)=1.0; end % state-to-output matrix in lagrange term
    Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
    W = eye(ny); % weight matrix in lagrange term
    for ii=1:nu W(ii,ii)=1e-2; end
    for ii=nu+1:nu+nx/2 W(ii,ii)=1e3; end
    for ii=nu+nx/2+1:nu+nx W(ii,ii)=1e-2; end
    W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
    yr = zeros(ny, 1); % output reference in lagrange term
    yr_e = zeros(ny_e, 1); % output reference in mayer term

    ocp_model.set('cost_Vu', Vu);
    ocp_model.set('cost_Vx', Vx);
    ocp_model.set('cost_Vx_e', Vx_e);
    ocp_model.set('cost_W', W);
    ocp_model.set('cost_W_e', W_e);
    ocp_model.set('cost_y_ref', yr);
    ocp_model.set('cost_y_ref_e', yr_e);
    ocp_model.set('dim_ny', ny);
    ocp_model.set('dim_ny_e', ny_e);
else % external, auto
    ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
    ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
end

% constraints
%Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,ii)=1.0; end
%lbx = -4*ones(nbx, 1);
%ubx =  4*ones(nbx, 1);
Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);


%% acados ocp model
ocp_model.set('name', model_name);
% dims
ocp_model.set('T', T);
ocp_model.set('dim_nx', nx);
ocp_model.set('dim_nu', nu);
ocp_model.set('dim_nbx', nbx);
ocp_model.set('dim_nbu', nbu);
% ocp_model.set('dim_ng', ng); % for general linear constraints
% ocp_model.set('dim_ng_e', ng_e);
ocp_model.set('dim_nh', nh);
ocp_model.set('dim_nh_e', nh_e);
% symbolics
ocp_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    ocp_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_xdot')
    ocp_model.set('sym_xdot', model.sym_xdot);
end

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.expr_f_expl);
else % irk irk_gnsf
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.expr_f_impl);
end
% constraints
ocp_model.set('constr_expr_h', model.expr_h);
ocp_model.set('constr_lh', lbu);
ocp_model.set('constr_uh', ubu);

x0 = [0; pi; 0; 0];
ocp_model.set('constr_x0', x0);

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
if (strcmp(nlp_solver, 'sqp')) % not available for sqp_rti
    ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
    ocp_opts.set('nlp_solver_tol_stat', tol);
    ocp_opts.set('nlp_solver_tol_eq', tol);
    ocp_opts.set('nlp_solver_tol_ineq', tol);
    ocp_opts.set('nlp_solver_tol_comp', tol);
end
ocp_opts.set('qp_solver', qp_solver);
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
ocp_opts.set('qp_solver_iter_max', qp_solver_iter_max);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', sim_method_num_steps);

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

% x_traj_init = zeros(nx, N+1);
x_traj_init = [linspace(0, 0, N+1); linspace(pi, 0, N+1); linspace(0, 0, N+1); linspace(0, 0, N+1)];
u_traj_init = zeros(nu, N);

%% prepare evaluation
n_executions = 100;
time_tot = zeros(n_executions,1);
time_lin = zeros(n_executions,1);
time_reg = zeros(n_executions,1);
time_qp_sol = zeros(n_executions,1);


%% call ocp solver in loop
for i=1:n_executions
    
    % initial state
    ocp.set('constr_x0', x0);

    % set trajectory initialization
    ocp.set('init_x', x_traj_init);
    ocp.set('init_u', u_traj_init);
    ocp.set('init_pi', zeros(nx, N))

    % solve
    ocp.solve();
    % get solution
    utraj = ocp.get('u');
    xtraj = ocp.get('x');

    %% evaluation
    status = ocp.get('status');
    sqp_iter = ocp.get('sqp_iter');
    time_tot(i) = ocp.get('time_tot');
    time_lin(i) = ocp.get('time_lin');
    time_reg(i) = ocp.get('time_reg');
    time_qp_sol(i) = ocp.get('time_qp_sol');

    if i == 1 || i == n_executions
        ocp.print('stat')
    end
end


%% plot average compuation times
if ~is_octave()
    time_total = sum(time_tot);
    time_linearize = sum(time_lin);
    time_regulariz = sum(time_reg);
    time_qp_solution = sum(time_qp_sol);

    figure;

    bar_vals = 1000 * [time_linearize; time_regulariz; time_qp_solution; ...
        time_total - time_linearize - time_regulariz - time_qp_solution] / n_executions;
    bar([1; nan], [bar_vals, nan(size(bar_vals))]' ,'stacked')
    legend('linearization', 'regularization', 'qp solution', 'remaining')
    ylabel('time in [ms]')
    title( [ strrep(cost_type, '_',' '), ' , sim: ' strrep(sim_method, '_',' '), ...
       ';  ', strrep(qp_solver, '_', ' ')] )
end

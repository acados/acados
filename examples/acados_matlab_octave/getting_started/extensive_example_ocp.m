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

addpath('../pendulum_on_cart_model')

check_acados_requirements()

print_level = 1;

%% discretization
N = 40;
T = 2.0; % time horizon length
h = T/N;

% nonuniform time grid
% N1 = 5;
% N2 = N - N1;
% time_steps = [( 1 * ones(N1,1)); 3 * ones(N2,1)];
% time_steps = T/sum(time_steps) * time_steps;

% uniform time grid
time_steps = T/N * ones(N,1);

shooting_nodes = zeros(N+1, 1);
for i = 1:N
    shooting_nodes(i+1) = sum(time_steps(1:i));
end

nlp_solver = 'sqp'; % sqp, sqp_rti
nlp_solver_exact_hessian = 'false';
regularize_method = 'convexify';
     % no_regularize, project, project_reduc_hess, mirror, convexify
nlp_solver_max_iter = 50;
tol = 1e-8;
qp_solver = 'partial_condensing_hpipm';
    % full_condensing_hpipm, partial_condensing_hpipm
    % full_condensing_qpoases, partial_condensing_osqp
qp_solver_cond_N = 5; % for partial condensing
qp_solver_cond_ric_alg = 0;
qp_solver_ric_alg = 0;
qp_solver_warm_start = 1; % 0: cold, 1: warm, 2: hot
qp_solver_iter_max = 1000; % default is 50; OSQP needs a lot sometimes.

% can vary for integrators
sim_method_num_stages = 1 * ones(N,1);
sim_method_num_steps = ones(N,1);
sim_method_num_stages(3:end) = 2;

%% model dynamics
model = pendulum_on_cart_model();

%% model to create the solver
ocp_model = acados_ocp_model();

%% dimensions
nx = model.nx;
nu = model.nu;

model_name = 'pendulum';

%% cost formulation
cost_formulation = 1;
switch cost_formulation
    case 1
        cost_type = 'linear_ls';
    case 2
        cost_type = 'ext_cost';
    otherwise
        cost_type = 'auto';
end

%% integrator type
integrator = 1;
switch integrator
    case 1
        sim_method = 'erk';
    case 2
        sim_method = 'irk';
    case 3
        sim_method = 'discrete';
    otherwise
        sim_method = 'irk_gnsf';
end

%% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);
if strcmp( cost_type, 'linear_ls' )
    ny = nu+nx; % number of outputs in lagrange term
    % input-to-output matrix in lagrange term
    Vu = zeros(ny, nu);
    Vu(1:nu,:) = eye(nu);
    % state-to-output matrix in lagrange term
    Vx = zeros(ny, nx);
    Vx(nu+1:end, :) = eye(nx);
    W = diag([1e-2, 1e3, 1e3, 1e-2, 1e-2]);

    % terminal cost term
    ny_e = nx; % number of outputs in terminal cost term
    Vx_e = eye(ny_e, nx);
    W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
    y_ref = zeros(ny, 1); % output reference in lagrange term
    y_ref_e = zeros(ny_e, 1); % output reference in mayer term

    ocp_model.set('cost_Vu', Vu);
    ocp_model.set('cost_Vx', Vx);
    ocp_model.set('cost_Vx_e', Vx_e);
    ocp_model.set('cost_W', W);
    ocp_model.set('cost_W_e', W_e);
    ocp_model.set('cost_y_ref', y_ref);
    ocp_model.set('cost_y_ref_e', y_ref_e);
else % external, auto
    ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
    ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
end

%% constraints
constraint_formulation_nonlinear = 0;
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);
if constraint_formulation_nonlinear % formulate constraint via h
    ocp_model.set('constr_expr_h', model.expr_h);
    ocp_model.set('constr_lh', lbu);
    ocp_model.set('constr_uh', ubu);
else % formulate constraint as bound on u
    Jbu = eye(nu);
    ocp_model.set('constr_Jbu', Jbu);
    ocp_model.set('constr_lbu', lbu);
    ocp_model.set('constr_ubu', ubu);
end

%% acados ocp model
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
if isfield(model, 'sym_z') % algebraic variables
    ocp_model.set('sym_z', model.sym_z);
end
if isfield(model, 'sym_p') % parameters
    ocp_model.set('sym_p', model.sym_p);
end

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.expr_f_expl);
elseif (strcmp(sim_method, 'irk') || strcmp(sim_method, 'irk_gnsf'))
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.expr_f_impl);
elseif strcmp(sim_method, 'discrete')
    ocp_model.set('dyn_type', 'discrete');
    % build explicit euler discrete integrator
    import casadi.*
    expl_ode_fun = Function([model_name,'_expl_ode_fun'], ...
            {model.sym_x, model.sym_u}, {model.expr_f_expl});
    dyn_expr_phi = model.sym_x + T/N * expl_ode_fun(model.sym_x, model.sym_u);
    ocp_model.set('dyn_expr_phi', dyn_expr_phi)
    if ~all(time_steps == T/N)
        disp('nonuniform time discretization with discrete dynamics should not be used');
        keyboard
    end
end

x0 = [0; pi; 0; 0];
ocp_model.set('constr_x0', x0);

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
if (exist('time_steps', 'var'))
	ocp_opts.set('time_steps', time_steps);
end

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

ocp_opts.set('exact_hess_dyn', 1);
ocp_opts.set('exact_hess_cost', 1);
ocp_opts.set('exact_hess_constr', 1);

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

x_traj_init = zeros(nx, N+1);
x_traj_init(2, :) = linspace(pi, 0, N+1); % initialize theta

u_traj_init = zeros(nu, N);

%% prepare evaluation
n_executions = 1;
time_tot = zeros(n_executions,1);
time_lin = zeros(n_executions,1);
time_reg = zeros(n_executions,1);
time_qp_sol = zeros(n_executions,1);

ocp.set('print_level', print_level)

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

% get slack values
for i = 0:N-1
    sl = ocp.get('sl', i);
    su = ocp.get('su', i);
    t = ocp.get('t', i);
end
sl = ocp.get('sl', N);
su = ocp.get('su', N);

% get cost value
cost_val_ocp = ocp.get_cost();


%% get QP matrices:
% See https://docs.acados.org/problem_formulation
% either stage wise
for stage = [0, N-1]
    % Note loop over field doesnt work because stupid matlab diff between
    % chars and strings
    % dynamics
    field = 'qp_A';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_B';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_b';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)

    % cost
    field = 'qp_R';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_Q';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_r';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)

    % constraints
    field = 'qp_C';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_D';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)

    field = 'qp_lg';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_ug';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)

    field = 'qp_lbx';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_ubx';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)

    field = 'qp_lbu';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
    field = 'qp_ubu';
    disp(strcat(field, " at stage ", num2str(stage), " = "));
    ocp.get(field, stage)
end

stage = N;
field = 'qp_Q';
disp(strcat(field, " at stage ", num2str(stage), " = "));
ocp.get(field, stage)
field = 'qp_R';
disp(strcat(field, " at stage ", num2str(stage), " = "));
ocp.get(field, stage)

% or for all stages
qp_Q = ocp.get('qp_Q');

%% Plot trajectories
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
for i=1:length(States)
    subplot(length(States), 1, i);
    plot(shooting_nodes, xtraj(i,:)); grid on;
    ylabel(States{i});
    xlabel('t [s]')
end

figure
stairs(shooting_nodes, [utraj'; utraj(end)])

ylabel('F [N]')
xlabel('t [s]')
grid on
if is_octave()
    waitforbuttonpress;
end

%% plot average compuation times
% if ~is_octave()
%     time_total = sum(time_tot);
%     time_linearize = sum(time_lin);
%     time_regulariz = sum(time_reg);
%     time_qp_solution = sum(time_qp_sol);
% 
%     figure;
% 
%     bar_vals = 1000 * [time_linearize; time_regulariz; time_qp_solution; ...
%         time_total - time_linearize - time_regulariz - time_qp_solution] / n_executions;
%     bar([1; nan], [bar_vals, nan(size(bar_vals))]' ,'stacked')
%     legend('linearization', 'regularization', 'qp solution', 'remaining')
%     ylabel('time in [ms]')
%     title( [ strrep(cost_type, '_',' '), ' , sim: ' strrep(sim_method, '_',' '), ...
%        ';  ', strrep(qp_solver, '_', ' ')] )
% end

%% test templated solver
% if ~ispc()
%     % MEX wrapper around templated solver not supported for Windows yet
%     % it can be used and tested via Simulink though.
%     disp('testing templated solver');
%     ocp.generate_c_code;
%     cd c_generated_code/
%     command = strcat('t_ocp = ', model_name, '_mex_solver');
%     eval( command );
%
%     t_ocp.set('print_level', print_level)
%
%     % initial state
%     t_ocp.set('constr_x0', x0);
%
%     % set trajectory initialization
%     t_ocp.set('init_x', x_traj_init);
%     t_ocp.set('init_u', u_traj_init);
%     t_ocp.set('init_pi', zeros(nx, N))
%
%     t_ocp.solve()
%     xt_traj = t_ocp.get('x');
%     ut_traj = t_ocp.get('u');
%
%     error_X_mex_vs_mex_template = max(max(abs(xt_traj - xtraj)))
%     error_U_mex_vs_mex_template = max(max(abs(ut_traj - utraj)))
%
%     t_ocp.print('stat')
%     cost_val_t_ocp = t_ocp.get_cost();
%     clear t_ocp
%     cd ..
% end

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

addpath('../wind_turbine_nx6/');

% qpOASES too slow to test on CI
for itest = 1:3

%% arguments
compile_interface = 'auto';
codgen_model = 'true';
% simulation
sim_method = 'irk';
sim_sens_forw = 'false';
sim_num_stages = 4;
sim_num_steps = 1;
% ocp
ocp_N = 40;
ocp_nlp_solver = 'sqp';
%ocp_nlp_solver = 'sqp_rti';
ocp_nlp_solver_exact_hessian = 'false';
%ocp_nlp_solver_exact_hessian = 'true';
regularize_method = 'no_regularize';
%regularize_method = 'project';
%regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
ocp_nlp_solver_max_iter = 50;
ocp_nlp_solver_tol_stat = 1e-8;
ocp_nlp_solver_tol_eq   = 1e-8;
ocp_nlp_solver_tol_ineq = 1e-8;
ocp_nlp_solver_tol_comp = 1e-8;
ocp_nlp_solver_ext_qp_res = 1;
switch itest
case 1
    ocp_qp_solver = 'partial_condensing_hpipm';
case 2
    ocp_qp_solver = 'full_condensing_hpipm';
case 3
    ocp_qp_solver = 'full_condensing_daqp';
case 4
    ocp_qp_solver = 'full_condensing_qpoases';
end
fprintf(['\n\nrunning with qp solver ', ocp_qp_solver, '\n'])
ocp_qp_solver_cond_N = 5;
ocp_qp_solver_cond_ric_alg = 0;
ocp_qp_solver_ric_alg = 0;
ocp_qp_solver_warm_start = 2;
%ocp_sim_method = 'erk';
ocp_sim_method = 'irk';
ocp_sim_method_num_stages = 4 * ones(ocp_N, 1); % scalar or vector of size ocp_N;
ocp_sim_method_num_steps = 1 * ones(ocp_N, 1); % scalar or vector of size ocp_N;
ocp_sim_method_newton_iter = 3; % * ones(ocp_N, 1); % scalar or vector of size ocp_N;
%cost_type = 'linear_ls';
cost_type = 'nonlinear_ls';

% get references
compute_setup;


%% create model entries
model = ocp_model_wind_turbine_nx6;



%% dims
Ts = 0.2; % samplig time
T = ocp_N*Ts; %8.0; % horizon length time [s]
nx = model.nx; % 8
nu = model.nu; % 2
ny = 4; % number of outputs in lagrange term
ny_e = 2; % number of outputs in mayer term
nbx = 3;
nbu = nu;
nh = 1;
nh_e = 1;
ns = 2;
ns_e = 2;
%ns_e = 1;
nsbx = 1;
%nsbx = 0;
nsh = 1;
nsh_e = 1;
np = model.np; % 1

%% cost
% state-to-output matrix in lagrange term
Vx = zeros(ny, nx);
Vx(1, 1) = 1.0;
Vx(2, 5) = 1.0;
% input-to-output matrix in lagrange term
Vu = zeros(ny, nu);
Vu(3, 1) = 1.0;
Vu(4, 2) = 1.0;
% state-to-output matrix in mayer term
Vx_e = zeros(ny_e, nx);
Vx_e(1, 1) = 1.0;
Vx_e(2, 5) = 1.0;
% weight matrix in lagrange term
W = zeros(ny, ny);
W(1, 1) =  1.5114;
W(2, 1) = -0.0649;
W(1, 2) = -0.0649;
W(2, 2) =  0.0180;
W(3, 3) =  0.01;
W(4, 4) =  0.001;
% weight matrix in mayer term
W_e = zeros(ny_e, ny_e); 
W_e(1, 1) =  1.5114;
W_e(2, 1) = -0.0649;
W_e(1, 2) = -0.0649;
W_e(2, 2) =  0.0180;
% output reference in lagrange term
%yr = ... ;
% output reference in mayer term
%yr_e = ... ;
% slacks
Z = 1e2*eye(ns);
Z_e = 1e2*eye(ns_e);
z = 0e2*ones(ns,1);
z_e = 0e2*ones(ns_e,1);

%% constraints
% constants
dbeta_min = -8.0;
dbeta_max =  8.0;
dM_gen_min = -1.0;
dM_gen_max =  1.0;
OmegaR_min =  6.0/60*2*3.14159265359;
OmegaR_max = 13.0/60*2*3.14159265359;
beta_min =  0.0;
beta_max = 35.0;
M_gen_min = 0.0;
M_gen_max = 5.0;
Pel_min = 0.0;
Pel_max = 5.0; % 5.0

%acados_inf = 1e8;

% state bounds
Jbx = zeros(nbx, nx);
Jbx(1, 1) = 1.0;
Jbx(2, 7) = 1.0;
Jbx(3, 8) = 1.0;
lbx = [OmegaR_min; beta_min; M_gen_min];
ubx = [OmegaR_max; beta_max; M_gen_max];
% input bounds
Jbu = eye(nu);
lbu = [dbeta_min; dM_gen_min];
ubu = [dbeta_max; dM_gen_max];
% nonlinear constraints (power constraint)
lh = Pel_min;
uh = Pel_max;
lh_e = Pel_min;
uh_e = Pel_max;
% soft box state constraints
Jsbx = zeros(nbx, nsbx);
Jsbx(1, 1) = 1.0;
% soft nonlinear constraints
Jsh = zeros(nh, nsh);
Jsh(1, 1) = 1.0;
Jsh_e = zeros(nh_e, nsh_e);
Jsh_e(1, 1) = 1.0;


%% acados ocp model
ocp_model = acados_ocp_model();
ocp_model.set('T', T);

%% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);
ocp_model.set('sym_p', model.sym_p);
%% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);
if (strcmp(cost_type, 'linear_ls'))
    ocp_model.set('cost_Vu', Vu);
    ocp_model.set('cost_Vx', Vx);
    ocp_model.set('cost_Vx_e', Vx_e);
else % nonlinear_ls
    ocp_model.set('cost_expr_y', model.expr_y);
    ocp_model.set('cost_expr_y_e', model.expr_y_e);
end
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', W_e);
ocp_model.set('cost_Z', Z);
ocp_model.set('cost_Z_e', Z_e);
ocp_model.set('cost_z', z);
ocp_model.set('cost_z_e', z_e);
%% dynamics
if (strcmp(ocp_sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.expr_f_expl);
else % irk
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.expr_f_impl);
end
%% constraints
% state bounds
ocp_model.set('constr_Jbx', Jbx);
ocp_model.set('constr_lbx', lbx);
ocp_model.set('constr_ubx', ubx);
ocp_model.set('constr_Jbx_e', Jbx);
ocp_model.set('constr_lbx_e', lbx);
ocp_model.set('constr_ubx_e', ubx);
% input bounds
ocp_model.set('constr_Jbu', Jbu);
ocp_model.set('constr_lbu', lbu);
ocp_model.set('constr_ubu', ubu);
% nonlinear constraints
ocp_model.set('constr_expr_h', model.expr_h);
ocp_model.set('constr_lh', lh);
ocp_model.set('constr_uh', uh);
ocp_model.set('constr_expr_h_e', model.expr_h_e);
ocp_model.set('constr_lh_e', lh_e);
ocp_model.set('constr_uh_e', uh_e);
% soft nonlinear constraints
ocp_model.set('constr_Jsbx', Jsbx);
ocp_model.set('constr_Jsbx_e', Jsbx);
ocp_model.set('constr_Jsh', Jsh);
ocp_model.set('constr_Jsh_e', Jsh_e);

% initial state dummy
ocp_model.set('constr_x0', x0_ref);
%
ocp_model.model_struct;



%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_interface', compile_interface);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme_N', ocp_N);
ocp_opts.set('nlp_solver', ocp_nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', ocp_nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', ocp_nlp_solver_ext_qp_res);
if (strcmp(ocp_nlp_solver, 'sqp'))
    ocp_opts.set('nlp_solver_max_iter', ocp_nlp_solver_max_iter);
    ocp_opts.set('nlp_solver_tol_stat', ocp_nlp_solver_tol_stat);
    ocp_opts.set('nlp_solver_tol_eq', ocp_nlp_solver_tol_eq);
    ocp_opts.set('nlp_solver_tol_ineq', ocp_nlp_solver_tol_ineq);
    ocp_opts.set('nlp_solver_tol_comp', ocp_nlp_solver_tol_comp);
end
ocp_opts.set('qp_solver', ocp_qp_solver);
ocp_opts.set('qp_solver_iter_max', 500);
if (strcmp(ocp_qp_solver, 'partial_condensing_hpipm'))
    ocp_opts.set('qp_solver_cond_N', ocp_qp_solver_cond_N);
    ocp_opts.set('qp_solver_cond_ric_alg', ocp_qp_solver_cond_ric_alg);
    ocp_opts.set('qp_solver_ric_alg', ocp_qp_solver_ric_alg);
    ocp_opts.set('qp_solver_warm_start', ocp_qp_solver_warm_start);
end
ocp_opts.set('sim_method', ocp_sim_method);
ocp_opts.set('sim_method_num_stages', ocp_sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', ocp_sim_method_num_steps);
ocp_opts.set('sim_method_newton_iter', ocp_sim_method_newton_iter);
ocp_opts.set('regularize_method', 'no_regularize');
ocp_opts.set('ext_fun_compile_flags', '');

ocp_opts.set('parameter_values', wind0_ref(:,1));

ocp_opts.opts_struct;



%% acados ocp
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);
%ocp
%ocp.C_ocp

%% acados sim model
sim_model = acados_sim_model();
% symbolics
sim_model.set('sym_x', model.sym_x);
if isfield(model, 'sym_u')
    sim_model.set('sym_u', model.sym_u);
end
if isfield(model, 'sym_xdot')
    sim_model.set('sym_xdot', model.sym_xdot);
end
if isfield(model, 'sym_p')
    sim_model.set('sym_p', model.sym_p);
end
% model
sim_model.set('T', T/ocp_N);
if (strcmp(sim_method, 'erk'))
    sim_model.set('dyn_type', 'explicit');
    sim_model.set('dyn_expr_f', model.expr_f_expl);
else % irk
    sim_model.set('dyn_type', 'implicit');
    sim_model.set('dyn_expr_f', model.expr_f_impl);
end



%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', sim_num_stages);
sim_opts.set('num_steps', sim_num_steps);
sim_opts.set('method', sim_method);
sim_opts.set('sens_forw', sim_sens_forw);

%sim_opts.opts_struct



%% acados sim
% create sim
sim = acados_sim(sim_model, sim_opts);


%% closed loop simulation
n_sim = 100;
n_sim_max = length(wind0_ref) - ocp_N;
if n_sim>n_sim_max
    n_sim = s_sim_max;
end
x_sim = zeros(nx, n_sim+1);
x_sim(:,1) = x0_ref; % initial state
u_sim = zeros(nu, n_sim);

sqp_iter_sim = zeros(n_sim,1);
time_ext = zeros(n_sim, 1);
time_tot = zeros(n_sim, 1);
time_lin = zeros(n_sim, 1);
time_qp_sol = zeros(n_sim, 1);

% set trajectory initialization
x_traj_init = repmat(x0_ref, 1, ocp_N+1);
u_traj_init = repmat(u0_ref, 1, ocp_N);
pi_traj_init = zeros(nx, ocp_N);


for ii=1:n_sim

%    fprintf('\nsimulation step %d\n', ii);

    tic

    % set x0
    ocp.set('constr_x0', x_sim(:,ii));
    % set parameter
    for jj=0:ocp_N-1
        ocp.set('p', wind0_ref(:,ii+jj), jj);
    end

    % set reference (different at each stage)
    for jj=0:ocp_N-1
        ocp.set('cost_y_ref', y_ref(:,ii+jj), jj);
    end
    ocp.set('cost_y_ref_e', y_ref(1:ny_e,ii+ocp_N));

    % set trajectory initialization (if not, set internally using previous solution)
    ocp.set('init_x', x_traj_init);
    ocp.set('init_u', u_traj_init);
    ocp.set('init_pi', pi_traj_init);

    % solve
    ocp.solve();

    % get solution
    x = ocp.get('x');
    u = ocp.get('u');
    pi = ocp.get('pi');

    % store first input
    u_sim(:,ii) = ocp.get('u', 0);

    % set initial state of sim
    sim.set('x', x_sim(:,ii));
    % set input in sim
    sim.set('u', u_sim(:,ii));
    % set parameter
    sim.set('p', wind0_ref(:,ii));

    % simulate state
    sim.solve();

    % get new state
    x_sim(:,ii+1) = sim.get('xn');

    % shift trajectory for initialization
%    x_traj_init = [x(:,2:ocp_N+1), zeros(nx, 1)];
    x_traj_init = [x(:,2:ocp_N+1), x(:,ocp_N+1)];
%    x_traj_init = [x(:,2:ocp_N+1), sim.get('xn')];
%    u_traj_init = [u(:,2:ocp_N), zeros(nu, 1)];
    u_traj_init = [u(:,2:ocp_N), u(:,ocp_N)];
    pi_traj_init = [pi(:,2:ocp_N), pi(:,ocp_N)];

    time_ext(ii) = toc;

    electrical_power = 0.944*97/100*x(1,1)*x(6,1);

    status = ocp.get('status');
    sqp_iter = ocp.get('sqp_iter');
    time_tot(ii) = ocp.get('time_tot');
    time_lin(ii) = ocp.get('time_lin');
    time_qp_sol(ii) = ocp.get('time_qp_sol');
    sqp_stats = ocp.get('stat');
    qp_iter = sqp_stats(:,7);

    sqp_iter_sim(ii) = sqp_iter;
    if status ~= 0
        ocp.print()
        error(['ocp_nlp solver returned status ', num2str(status), '!= 0 in simulation instance ', num2str(ii)]);
    end

    fprintf('\nstatus = %d, sqp_iter = %d, qp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms]), Pel = %f',...
            status, sqp_iter, sum(qp_iter), time_ext(ii)*1e3, time_tot(ii)*1e3, time_lin(ii)*1e3, time_qp_sol(ii)*1e3, electrical_power);

    if 0
        ocp.print('stat')
    end

end

% test setter
ocp.set('cost_z', ones(2,1), 1)
ocp.set('cost_Z', ones(2,1), 1)
ocp.set('cost_zl', ones(2,1), N-1)

% get slack values
for i = 0:N-1
    sl = ocp.get('sl', i);
    su = ocp.get('su', i);
    % test setters
    ocp.set('sl', sl, i);
    ocp.set('su', su, i);
    t = ocp.get('t', i);
end
sl = ocp.get('sl', N);
su = ocp.get('su', N);



electrical_power = 0.944*97/100*x_sim(1,:).*x_sim(6,:);

x_sim_ref = [   1.263425730522397
   0.007562725557589
  76.028289356099236
   0.007188510774546
   6.949049234224142
   3.892712459979240
   6.302629591941585
   3.882220255648666];

err_vs_ref = x_sim_ref - x_sim(:,end);

fprintf('\nmedian computation times: time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])', ...
        median(time_ext)*1e3, median(time_tot)*1e3, median(time_lin)*1e3, median(time_qp_sol)*1e3)

if status~=0
    error('test_ocp_wtnx6: solution failed!');
elseif err_vs_ref > 1e-14
    error('test_ocp_wtnx6: to high deviation from known result!');
elseif sqp_iter > 2
    error('test_ocp_wtnx6: sqp_iter > 2, this problem is typically solved within less iterations!');
else
    fprintf('\ntest_ocp_wtnx6: success!\n');
end

% figures
if 0
    figure;
    subplot(3,1,1);
    plot(0:n_sim, x_sim);
    xlim([0 n_sim]);
    ylabel('states');
    %legend('p', 'theta', 'v', 'omega');
    subplot(3,1,2);
    plot(0:n_sim-1, u_sim);
    xlim([0 n_sim]);
    ylabel('inputs');
    %legend('F');
    subplot(3,1,3);
    plot(0:n_sim, electrical_power);
    hold on
    plot([0 n_sim], [Pel_max Pel_max]);
    hold off
    xlim([0 n_sim]);
    ylim([4.0 6.0]);
    ylabel('electrical power');
    %legend('F');

    figure;
    plot(1:n_sim, sqp_iter_sim, 'rx');
    hold on
    plot([1 n_sim], [ocp_nlp_solver_max_iter ocp_nlp_solver_max_iter]);
    hold off
    ylim([0 ocp_nlp_solver_max_iter+1])
    ylabel('sqp iterations')
    xlabel('sqp calls')
    if is_octave()
        waitforbuttonpress;
    end
end

end

% remove temporary created files
delete('y_ref.mat')
delete('y_e_ref.mat')
delete('wind0_ref.mat')
delete('windN_ref.mat')

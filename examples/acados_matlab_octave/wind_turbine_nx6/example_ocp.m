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
% Wind turbine OCP (nx = 6/8) ported to the NEW acados MATLAB/Octave
% interface (>= v0.4.0), following the style of the linear_mass_spring
% reference example.

clear all
addpath('.');

% (optional) check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    warning('env.sh does not seem to be sourced. If the build fails, run: source env.sh');
end

%% 1. Load model (reuse the existing symbolic definition)
raw = ocp_model_wind_turbine_nx6();

nx = raw.nx;   % 8
nu = raw.nu;   % 2
np = raw.np;   % 1

model = AcadosModel();
model.name        = 'wind_turbine_nx6';
model.x           = raw.sym_x;
model.u           = raw.sym_u;
model.xdot        = raw.sym_xdot;
model.p           = raw.sym_p;
model.f_expl_expr = raw.expr_f_expl;
model.f_impl_expr = raw.expr_f_impl;

%% 2. Set up OCP object
ocp = AcadosOcp();
ocp.model = model;

% --- Horizon & shooting nodes ---
N  = 40;
Ts = 0.2;        % sampling time
T  = N * Ts;     % horizon length [s]

%% 3. Cost  (linear least squares)
ny   = 4;   % number of outputs in Lagrange term
ny_e = 2;   % number of outputs in Mayer term

ocp.cost.cost_type   = 'LINEAR_LS';
ocp.cost.cost_type_e = 'LINEAR_LS';
ocp.cost.cost_type_0 = 'LINEAR_LS';

% state-to-output matrix (Lagrange)
Vx = zeros(ny, nx);
Vx(1, 1) = 1.0;
Vx(2, 5) = 1.0;
% input-to-output matrix (Lagrange)
Vu = zeros(ny, nu);
Vu(3, 1) = 1.0;
Vu(4, 2) = 1.0;
% state-to-output matrix (Mayer)
Vx_e = zeros(ny_e, nx);
Vx_e(1, 1) = 1.0;
Vx_e(2, 5) = 1.0;

% weight matrix (Lagrange)
W = zeros(ny, ny);
W(1, 1) =  1.5114;
W(2, 1) = -0.0649;
W(1, 2) = -0.0649;
W(2, 2) =  0.0180;
W(3, 3) =  0.01;
W(4, 4) =  0.001;
% weight matrix (Mayer)
W_e = zeros(ny_e, ny_e);
W_e(1, 1) =  1.5114;
W_e(2, 1) = -0.0649;
W_e(1, 2) = -0.0649;
W_e(2, 2) =  0.0180;

% initial-node cost = same as intermediate (cost_type_0 = LINEAR_LS)
ocp.cost.Vx_0   = Vx;
ocp.cost.Vu_0   = Vu;
ocp.cost.W_0    = W;
ocp.cost.yref_0 = zeros(ny, 1);   % overwritten at runtime

ocp.cost.Vx     = Vx;
ocp.cost.Vu     = Vu;
ocp.cost.W      = W;
ocp.cost.yref   = zeros(ny, 1);   % overwritten at runtime

ocp.cost.Vx_e   = Vx_e;
ocp.cost.W_e    = W_e;
ocp.cost.yref_e = zeros(ny_e, 1); % overwritten at runtime

% --- slack penalties ---
Z   = 1e2;  z   = 0e2;
Z_e = 1e2;  z_e = 0e2;

% intermediate soft nonlinear constraint
ocp.cost.Zl = Z;  ocp.cost.Zu = Z;
ocp.cost.zl = z;  ocp.cost.zu = z;
% terminal soft nonlinear constraint
ocp.cost.Zl_e = Z_e;  ocp.cost.Zu_e = Z_e;
ocp.cost.zl_e = z_e;  ocp.cost.zu_e = z_e;

use_soft_h0_constraint = 1;
if use_soft_h0_constraint
    % initial-node soft nonlinear constraint
    ocp.cost.Zl_0 = Z;  ocp.cost.Zu_0 = Z;
    ocp.cost.zl_0 = z;  ocp.cost.zu_0 = z;
end

%% 4. Constraints
% --- bound constants ---
dbeta_min  = -8.0;   dbeta_max  =  8.0;
dM_gen_min = -1.0;   dM_gen_max =  1.0;
OmegaR_min =  6.0/60*2*pi;   OmegaR_max = 13.0/60*2*pi;
beta_min   =  0.0;   beta_max   = 35.0;
M_gen_min  =  0.0;   M_gen_max  =  5.0;
Pel_min    =  0.0;   Pel_max    =  5.0;

% --- state bounds (Jbx -> idxbx, 0-based indices) ---
% old Jbx selected x(1), x(7), x(8)  ->  idxbx = [0, 6, 7]
ocp.constraints.idxbx = [0; 6; 7];
ocp.constraints.lbx   = [OmegaR_min; beta_min; M_gen_min];
ocp.constraints.ubx   = [OmegaR_max; beta_max; M_gen_max];

% --- input bounds (Jbu = eye(nu) -> idxbu = [0, 1]) ---
ocp.constraints.idxbu = [0; 1];
ocp.constraints.lbu   = [dbeta_min; dM_gen_min];
ocp.constraints.ubu   = [dbeta_max; dM_gen_max];

% --- nonlinear power constraint h, on initial / intermediate / terminal ---
ocp.model.con_h_expr   = raw.expr_h;
ocp.constraints.lh     = Pel_min;
ocp.constraints.uh     = Pel_max;

ocp.model.con_h_expr_e = raw.expr_h_e;
ocp.constraints.lh_e   = Pel_min;
ocp.constraints.uh_e   = Pel_max;

if use_soft_h0_constraint
    ocp.model.con_h_expr_0 = raw.expr_h;
    ocp.constraints.lh_0   = Pel_min;
    ocp.constraints.uh_0   = Pel_max;
end

% --- soft nonlinear constraints (Jsh -> idxsh, 0-based) ---
% lsh/ush (slack lower/upper bounds) default to 0 and are filled in
ocp.constraints.idxsh   = 0;   % soften the (only) h constraint
ocp.constraints.idxsh_e = 0;
if use_soft_h0_constraint
    ocp.constraints.idxsh_0 = 0;
end

% --- (dummy) initial state, overwritten at runtime ---
ocp.constraints.x0 = zeros(nx, 1);

% --- parameter default value (np = 1), overwritten per stage at runtime ---
ocp.parameter_values = zeros(np, 1);

%% 5. Create solver
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf        = T;

% --- Solver options (mapping from the old ocp_opts) ---
ocp.solver_options.nlp_solver_type       = 'SQP';            % was 'sqp'
ocp.solver_options.hessian_approx        = 'GAUSS_NEWTON';   % nlp_solver_exact_hessian = 'false'
ocp.solver_options.regularize_method     = 'NO_REGULARIZE';  % was 'no_regularize'
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.nlp_solver_max_iter   = 200;

ocp.solver_options.qp_solver             = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N      = 5;
ocp.solver_options.qp_solver_cond_ric_alg = 0;
ocp.solver_options.qp_solver_ric_alg      = 0;
ocp.solver_options.qp_solver_warm_start   = 0;
ocp.solver_options.qp_solver_iter_max     = 50;   % old: qp_solver_max_iter

% --- Integrator (was sim_method = 'irk') ---
ocp.solver_options.integrator_type    = 'IRK';    % use 'ERK' with f_expl_expr instead
ocp.solver_options.sim_method_num_stages = 4;
ocp.solver_options.sim_method_num_steps  = 1;     % (was 4 for erk, 1 for irk)

ocp_solver = AcadosOcpSolver(ocp);

%% 6. References / initialization
% provides x0_ref, u0_ref, wind0_ref, y_ref
compute_setup;

% trajectory initialization
x_traj_init = repmat(x0_ref, 1, N+1);
u_traj_init = repmat(u0_ref, 1, N);
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% set x0
ocp_solver.set('constr_x0', x0_ref);

% set time-varying parameter (wind) per stage
for jj = 0:N
    ocp_solver.set('p', wind0_ref(:, jj+1), jj);
end

% set references per stage
for jj = 0:N-1
    ocp_solver.set('cost_y_ref', y_ref(:, jj+1), jj);
end
ocp_solver.set('cost_y_ref', y_ref(1:ny_e, N+1), N);

%% 7. Solve
disp('before solve')
tic;
ocp_solver.solve();
time_ext = toc;
disp('after solve')

% get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');

electrical_power = 0.944*97/100 * x(1,:) .* x(6,:);

%% 8. Diagnostics
status      = ocp_solver.get('status');
sqp_iter    = ocp_solver.get('sqp_iter');
time_tot    = ocp_solver.get('time_tot');
time_lin    = ocp_solver.get('time_lin');
time_reg    = ocp_solver.get('time_reg');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf(['\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] ', ...
         '(time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n'], ...
        status, sqp_iter, time_ext*1e3, time_tot*1e3, ...
        time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

ocp_solver.print('stat');

%% 9. Figures
figure;
subplot(3,1,1);
plot(0:N, x);
xlim([0 N]); ylabel('states');
subplot(3,1,2);
plot(0:N-1, u);
xlim([0 N]); ylabel('inputs');
subplot(3,1,3);
plot(0:N, electrical_power);
hold on
plot([0 N], [Pel_max Pel_max]);
hold off
xlim([0 N]); ylim([4.5 6.0]);
ylabel('electrical power');

stat = ocp_solver.get('stat');
figure;
plot(0:size(stat,1)-1, log10(stat(:,2)), 'r-x'); hold on
plot(0:size(stat,1)-1, log10(stat(:,3)), 'b-x');
plot(0:size(stat,1)-1, log10(stat(:,4)), 'g-x');
plot(0:size(stat,1)-1, log10(stat(:,5)), 'k-x'); hold off
xlabel('iter'); ylabel('res')

if status == 0
    fprintf('\nsuccess!\n\n');
else
    fprintf('\nsolution failed!\n\n');
end

if is_octave()
    waitforbuttonpress;
end
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

model_name = 'simple_dae';

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% options
compile_interface = 'auto'; % true, false

% ocp
N = 20;
nlp_solver = 'sqp'; % sqp, sqp_rti
nlp_solver_exact_hessian = 'true';
%nlp_solver_exact_hessian = 'true';
regularize_method = 'no_regularize';
%regularize_method = 'project_reduc_hess';
nlp_solver_max_iter = 10;
nlp_solver_tol_stat = 1e-12;
nlp_solver_tol_eq   = 1e-12;
nlp_solver_tol_ineq = 1e-12;
nlp_solver_tol_comp = 1e-12;
nlp_solver_ext_qp_res = 1;
nlp_solver_step_length = 0.7;
%qp_solver = 'full_condensing_qpoases'; % partial_condensing_hpipm
qp_solver = 'partial_condensing_hpipm'; % partial_condensing_hpipm
qp_solver_cond_N = 5;
qp_solver_warm_start = 0;
qp_solver_cond_ric_alg = 1; % 0: dont factorize hessian in the condensing; 1: factorize
qp_solver_ric_alg = 1; % HPIPM specific
%ocp_sim_method = 'irk'; % irk, irk_gnsf
ocp_sim_method = 'irk'; % irk, irk_gnsf
ocp_sim_method_num_stages = 6; % scalar or vector of size ocp_N;
ocp_sim_method_num_steps = 4; % scalar or vector of size ocp_N;
ocp_sim_method_newton_iter = 3; % scalar or vector of size ocp_N;

% selectors for example variants
constr_variant = 1; % 0: x bounds; 1: z bounds
cost_variant = 1; % 0: ls on u,x; 1: ls on u,z; (not implemented yet: 2: nls on u,z)

% get model
model = simple_dae_model();

nx = length(model.sym_x);
nu = length(model.sym_u);
nz = length(model.sym_z);
ny = nx+nu;
ny_e = nx;

T = 1.0;
h = T/N;

Wu = 1e-3*eye(nu);
Wx = 1e1*eye(nx);
W = [Wu, zeros(nu,nx); zeros(nx,nu), Wx];
Vu = [eye(nu); zeros(nx,nu)];
Vx = [zeros(nu,nx); eye(nx)];
Vx_e = eye(nx);

lb = [2; -2];
ub = [4;  2];

x0 = [3; -1.8];


%% acados ocp model
ocp_model = acados_ocp_model();

ocp_model.set('T', T);
ocp_model.set('name', model_name);

% symbolics
ocp_model.set('sym_x', model.sym_x);
ocp_model.set('sym_u', model.sym_u);
ocp_model.set('sym_xdot', model.sym_xdot);
ocp_model.set('sym_z', model.sym_z);

% cost
if cost_variant==0
    ocp_model.set('cost_type', 'linear_ls');
    ocp_model.set('cost_Vu', Vu);
    ocp_model.set('cost_Vx', Vx);
elseif cost_variant==1
    ocp_model.set('cost_type', 'linear_ls');
    ocp_model.set('cost_Vu', Vu);
    ocp_model.set('cost_Vz', Vx);
    ocp_model.set('cost_Vx', zeros(ny, nx));
else
    ocp_model.set('cost_type', 'nonlinear_ls');
    ocp_model.set('cost_expr_y', model.expr_y);
end
ocp_model.set('cost_type_e', 'linear_ls');
ocp_model.set('cost_Vx_e', Vx_e);
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', Wx);
%ocp_model.set('cost_y_ref', yr);
%ocp_model.set('cost_y_ref_e', yr_e);

% dynamics
ocp_model.set('dyn_type', 'implicit');
ocp_model.set('dyn_expr_f', model.expr_f_impl);

% constraints
ocp_model.set('constr_x0', x0);
if constr_variant==0
    ocp_model.set('constr_Jbx', eye(nx));
    ocp_model.set('constr_lbx', lb);
    ocp_model.set('constr_ubx', ub);
else
    ocp_model.set('constr_expr_h', model.expr_h);
    ocp_model.set('constr_lh', lb);
    ocp_model.set('constr_uh', ub);
    ocp_model.set('constr_expr_h_e', model.expr_h_e);
    ocp_model.set('constr_lh_e', lb);
    ocp_model.set('constr_uh_e', ub);
end



%% acados ocp opts
ocp_opts = acados_ocp_opts();

ocp_opts.set('compile_interface', compile_interface);
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);
ocp_opts.set('nlp_solver_step_length', nlp_solver_step_length);
if (strcmp(nlp_solver, 'sqp'))
    ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
    ocp_opts.set('nlp_solver_tol_stat', nlp_solver_tol_stat);
    ocp_opts.set('nlp_solver_tol_eq', nlp_solver_tol_eq);
    ocp_opts.set('nlp_solver_tol_ineq', nlp_solver_tol_ineq);
    ocp_opts.set('nlp_solver_tol_comp', nlp_solver_tol_comp);
end
ocp_opts.set('qp_solver', qp_solver);
% overwrite default qp solver tol which is same as nlp tol
%ocp_opts.set('qp_solver_tol_stat', qp_solver_tol_stat);
%ocp_opts.set('qp_solver_tol_eq', qp_solver_tol_eq);
%ocp_opts.set('qp_solver_tol_ineq', qp_solver_tol_ineq);
%ocp_opts.set('qp_solver_tol_comp', qp_solver_tol_comp);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
    ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
    ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
    ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
    ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
end
ocp_opts.set('sim_method', ocp_sim_method);
ocp_opts.set('sim_method_num_stages', ocp_sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', ocp_sim_method_num_steps);
ocp_opts.set('sim_method_newton_iter', ocp_sim_method_newton_iter);

ocp_opts.set('exact_hess_dyn', 1);
ocp_opts.set('exact_hess_cost', 1);
ocp_opts.set('exact_hess_constr', 1);



%% acados ocp
ocp_solver = acados_ocp(ocp_model, ocp_opts);

ocp_solver.solve();

stat = ocp_solver.get('stat');

ocp_solver.print('stat')

status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
sqp_time = ocp_solver.get('time_tot');
%if status ~= 0
%    keyboard
%end

format short e
% get solution, e.g. for initialization of next NLP
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
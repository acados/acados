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

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% arguments
compile_interface = 'auto'; %'auto';
codgen_model = 'true';
gnsf_detect_struct = 'true';

% discretization
N = 20;
T = 1; % horizon length time
h = T/N;

nlp_solver = 'sqp';
%nlp_solver = 'sqp_rti';
%nlp_solver_exact_hessian = 'false';
nlp_solver_exact_hessian = 'true';
%regularize_method = 'no_regularize';
%regularize_method = 'project';
regularize_method = 'project_reduc_hess';
%regularize_method = 'mirror';
%regularize_method = 'convexify';
nlp_solver_max_iter = 100; %100;
nlp_solver_tol_stat = 1e-8;
nlp_solver_tol_eq   = 1e-8;
nlp_solver_tol_ineq = 1e-8;
nlp_solver_tol_comp = 1e-8;
nlp_solver_ext_qp_res = 1;
qp_solver = 'partial_condensing_hpipm';
%qp_solver = 'full_condensing_hpipm';
%qp_solver = 'full_condensing_qpoases';
qp_solver_cond_N = 5;
qp_solver_cond_ric_alg = 0;
qp_solver_ric_alg = 0;
qp_solver_warm_start = 0;
qp_solver_max_iter = 100;
%sim_method = 'erk';
sim_method = 'irk';
%sim_method = 'irk_gnsf';
sim_method_num_stages = 4;
sim_method_num_steps = 3;
cost_type = 'linear_ls';
%cost_type = 'ext_cost';
model_name = 'ocp_pendulum';


%% create model entries
model = pendulum_on_cart_model();

% dims
nx = model.nx;
nu = model.nu;
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
if 0
    nbx = 0;
    nbu = nu;
    ng = 0;
    ng_e = 0;
    nh = 0;
    nh_e = 0;
else
    nbx = 0;
    nbu = 0;
    ng = 0;
    ng_e = 0;
    nh = nu;
    nh_e = 0;
end

% cost
% input-to-output matrix in lagrange term
Vu = zeros(ny, nu);
Vu(1:nu,:) = eye(nu);
% state-to-output matrix in lagrange term
Vx = zeros(ny, nx);
Vx(nu+1:end, :) = eye(nx);
% W = diag([1e-2, 1e3, 1e3, 1e-2, 1e-2]);
% high penalty on u -> no active constraints
W = diag([1e0, 1e3, 1e3, 1e-2, 1e-2]);

% terminal cost term
ny_e = nx; % number of outputs in terminal cost term
Vx_e = eye(ny_e, nx);
W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
y_ref = zeros(ny, 1); % output reference in lagrange term
y_ref_e = zeros(ny_e, 1); % output reference in mayer term

% constraints
x0 = [0; pi; 0; 0];
%Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,ii)=1.0; end
%lbx = -4*ones(nbx, 1);
%ubx =  4*ones(nbx, 1);
Jbu = eye(nbu, nu);
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);


%% acados ocp model
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
%if (strcmp(cost_type, 'linear_ls'))
    ocp_model.set('cost_Vu', Vu);
    ocp_model.set('cost_Vx', Vx);
    ocp_model.set('cost_Vx_e', Vx_e);
    ocp_model.set('cost_W', W);
    ocp_model.set('cost_W_e', W_e);
    ocp_model.set('cost_y_ref', y_ref);
    ocp_model.set('cost_y_ref_e', y_ref_e);
%else % if (strcmp(cost_type, 'ext_cost'))
%    ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
%    ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
%end

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', model.expr_f_expl);
else % irk irk_gnsf
    ocp_model.set('dyn_type', 'implicit');
    ocp_model.set('dyn_expr_f', model.expr_f_impl);
end
% constraints
ocp_model.set('constr_x0', x0);
if (nh>0)
    ocp_model.set('constr_expr_h', model.expr_h);
    ocp_model.set('constr_lh', lbu);
    ocp_model.set('constr_uh', ubu);
%    ocp_model.set('constr_expr_h_e', model.expr_h_e);
%    ocp_model.set('constr_lh_e', lh_e);
%    ocp_model.set('constr_uh_e', uh_e);
else
%    ocp_model.set('constr_Jbx', Jbx);
%    ocp_model.set('constr_lbx', lbx);
%    ocp_model.set('constr_ubx', ubx);
    ocp_model.set('constr_Jbu', Jbu);
    ocp_model.set('constr_lbu', lbu);
    ocp_model.set('constr_ubu', ubu);
end
% disp('ocp_model.model_struct')
% disp(ocp_model.model_struct)


%% acados ocp opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('compile_interface', compile_interface);
ocp_opts.set('codgen_model', codgen_model);
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
ocp_opts.set('regularize_method', regularize_method);
ocp_opts.set('nlp_solver_ext_qp_res', nlp_solver_ext_qp_res);
if (strcmp(nlp_solver, 'sqp'))
    ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
    ocp_opts.set('nlp_solver_tol_stat', nlp_solver_tol_stat);
    ocp_opts.set('nlp_solver_tol_eq', nlp_solver_tol_eq);
    ocp_opts.set('nlp_solver_tol_ineq', nlp_solver_tol_ineq);
    ocp_opts.set('nlp_solver_tol_comp', nlp_solver_tol_comp);
end
ocp_opts.set('qp_solver', qp_solver);
if (strcmp(qp_solver, 'partial_condensing_hpipm'))
    ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
    ocp_opts.set('qp_solver_ric_alg', qp_solver_ric_alg);
end
ocp_opts.set('qp_solver_cond_ric_alg', qp_solver_cond_ric_alg);
ocp_opts.set('qp_solver_warm_start', qp_solver_warm_start);
ocp_opts.set('qp_solver_iter_max', qp_solver_max_iter);
ocp_opts.set('sim_method', sim_method);
ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
ocp_opts.set('sim_method_num_steps', sim_method_num_steps);
if (strcmp(sim_method, 'irk_gnsf'))
    ocp_opts.set('gnsf_detect_struct', gnsf_detect_struct);
end

% disp('ocp_opts');
% disp(ocp_opts.opts_struct);


%% acados ocp
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);
% ocp
% disp('ocp.C_ocp');
% disp(ocp.C_ocp);
%ocp.model_struct

% set trajectory initialization
%x_traj_init = zeros(nx, N+1);
%for ii=1:N x_traj_init(:,ii) = [0; pi; 0; 0]; end
x_traj_init = [linspace(0, 0, N+1); linspace(pi, 0, N+1); linspace(0, 0, N+1); linspace(0, 0, N+1)];

u_traj_init = zeros(nu, N);

% if not set, the trajectory is initialized with the previous solution
ocp.set('init_x', x_traj_init);
ocp.set('init_u', u_traj_init);

% solve ocp
tic;
ocp.solve();
time_ext = toc;

% get solution
u = ocp.get('u');
x = ocp.get('x');

%% evaluation
status = ocp.get('status');
sqp_iter = ocp.get('sqp_iter');
time_tot = ocp.get('time_tot');
time_lin = ocp.get('time_lin');
time_reg = ocp.get('time_reg');
time_qp_sol = ocp.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n',...
    status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

ocp.print('stat');


%% figures
% plot trajectories
if 1
    for ii=1:N+1
        x_cur = x(:,ii);
    %    visualize;
    end

    figure;
    subplot(2,1,1);
    plot(0:N, x);
    title('trajectories')
    xlim([0 N]);
    legend('p', 'theta', 'v', 'omega');
    subplot(2,1,2);
    plot(0:N-1, u);
    xlim([0 N]);
    legend('F');
end

% plot residuals over iteraions
stat = ocp.get('stat');
if 0 && (strcmp(nlp_solver, 'sqp'))
    figure;
     plot([0: size(stat,1)-1], log10(stat(:,2)), 'r-x');
     hold on
     plot([0: size(stat,1)-1], log10(stat(:,3)), 'b-x');
     plot([0: size(stat,1)-1], log10(stat(:,4)), 'g-x');
     plot([0: size(stat,1)-1], log10(stat(:,5)), 'k-x');
%    semilogy(0: size(stat,1)-1, stat(:,2), 'r-x');
%    hold on
%    semilogy(0: size(stat,1)-1, stat(:,3), 'b-x');
%    semilogy(0: size(stat,1)-1, stat(:,4), 'g-x');
%    semilogy(0: size(stat,1)-1, stat(:,5), 'k-x');
    hold off
    xlabel('iter')
    ylabel('res')
    legend('res stat', 'res eq', 'res ineq', 'res compl');
end


if status==0
    fprintf('\nsuccess!\n\n');
else
    fprintf('\nsolution failed!\n\n');
end


%% paramteric sensitivity of solution
if 1
    field = 'ex'; % equality constraint on states
    stage = 0;
    index = 1;
    ocp.eval_param_sens(field, stage, index);

    sens_u = ocp.get('sens_u');
    sens_x = ocp.get('sens_x');

    % plot sensitivity
    figure
    subplot(2,1,1);
    plot(0:N, sens_x);
    title('sensitivities')
    xlim([0 N]);
    legend('p', 'theta', 'v', 'omega');
    subplot(2,1,2);
    plot(0:N-1, sens_u);
    xlim([0 N]);
    legend('F');

    % plot predicted solution
    figure
    subplot(2,1,1);
    plot(0:N, x+sens_x);
    title('predicted trajectories')
    xlim([0 N]);
    legend('p', 'theta', 'v', 'omega');
    subplot(2,1,2);
    plot(0:N-1, u+sens_u);
    xlim([0 N]);
    legend('F');

    for ii=1:N+1
        x_cur = x(:,ii)+sens_x(:,ii);
    %    visualize;
    end

end

sens_u = zeros(nx, N);
% get sensitivities w.r.t. initial state value with index
for index = 0:nx-1
    ocp.eval_param_sens(field, stage, index);
    sens_u(index+1,:) = ocp.get('sens_u');
end
disp('solution sensitivity dU_dx0')
disp(sens_u)


% qp_hess = ocp.get('qp_solver_cond_H');
% nv = size(qp_hess, 1);
% % make full
% for jj=1:nv
%     for ii=jj+1:nv
%         qp_hess(jj,ii) = qp_hess(ii,jj);
%     end
% end
% qp_hessian_cond_num = cond(qp_hess)


if is_octave()
    waitforbuttonpress;
end

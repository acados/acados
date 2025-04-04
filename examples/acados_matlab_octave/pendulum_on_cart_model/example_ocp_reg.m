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


clear all

check_acados_requirements()


%% model dynamics
model = get_pendulum_on_cart_model();
nx = length(model.x); % state size
nu = length(model.u); % input size

%% OCP formulation object
ocp = AcadosOcp();
ocp.model = model;

% discretization
N = 100;
h = 0.01;

nlp_solver_step_length = 1.0;
nlp_solver_exact_hessian = 'true';
regularize_method = 'PROJECT'; % PROJECT_REDUC_HESS, PROJECT, GERSHGORIN_LEVENBERG_MARQUARDT
nlp_solver_max_iter = 100; %100;
nlp_solver_tol_stat = 1e-8;
nlp_solver_tol_eq   = 1e-8;
nlp_solver_tol_ineq = 1e-8;
nlp_solver_tol_comp = 1e-8;

model_name = 'ocp_pendulum';

% dims
T = N*h; % horizon length time
ny = nu+nx; % number of outputs in path cost
ny_e = nx; % number of outputs in terminal cost term

% cost
Vu = zeros(ny, nu); for ii=1:nu Vu(ii,ii)=1.0; end % input-to-output matrix in lagrange term
Vx = zeros(ny, nx); for ii=1:nx Vx(nu+ii,ii)=1.0; end % state-to-output matrix in lagrange term
Vx_e = zeros(ny_e, nx); for ii=1:nx Vx_e(ii,ii)=1.0; end % state-to-output matrix in mayer term
W = eye(ny); % weight matrix in lagrange term
for ii=1:nu W(ii,ii)=1e-2; end
for ii=nu+1:nu+nx/2 W(ii,ii)=1e3; end
%for ii=nu+1:nu+nx/2 W(ii,ii)=1e1; end
for ii=nu+nx/2+1:nu+nx W(ii,ii)=1e-2; end
W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
yref = zeros(ny, 1); % output reference in lagrange term
yref_e = zeros(ny_e, 1); % output reference in mayer term

ocp.cost.cost_type = 'LINEAR_LS';
ocp.cost.Vx = Vx;
ocp.cost.Vu = Vu;
ocp.cost.Vx_e = Vx_e;
ocp.cost.W = W;
ocp.cost.W_e = W_e;
ocp.cost.yref = yref;
ocp.cost.yref_e = yref_e;

%
lbu = -80*ones(nu, 1);
ubu =  80*ones(nu, 1);

% constraints
x0 = [0; pi; 0; 0];
constr_expr_h = model.u;

% formulate as nonlinear constraint
ocp.constraints.x0 = x0;
ocp.model.con_h_expr_0 = constr_expr_h;
ocp.constraints.lh_0 = lbu;
ocp.constraints.uh_0 = ubu;
ocp.model.con_h_expr = constr_expr_h;
ocp.constraints.lh = lbu;
ocp.constraints.uh = ubu;

% options
ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = N;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.hessian_approx = 'EXACT';
ocp.solver_options.integrator_type = 'IRK';

ocp.solver_options.regularize_method = regularize_method;
ocp.solver_options.nlp_solver_ext_qp_res = 1;
ocp.solver_options.nlp_solver_step_length = nlp_solver_step_length;
ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter;
ocp.solver_options.nlp_solver_tol_stat = nlp_solver_tol_stat;
ocp.solver_options.nlp_solver_tol_eq = nlp_solver_tol_eq;
ocp.solver_options.nlp_solver_tol_ineq = nlp_solver_tol_ineq;
ocp.solver_options.nlp_solver_tol_comp = nlp_solver_tol_comp;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.qp_solver_cond_N = 5;
ocp.solver_options.qp_solver_ric_alg = 0;
ocp.solver_options.qp_solver_cond_ric_alg = 0;
ocp.solver_options.qp_solver_warm_start = 0;
ocp.solver_options.qp_solver_iter_max = 100;
ocp.solver_options.sim_method_num_stages = 2;
ocp.solver_options.sim_method_num_steps = 1;

%% create solver
ocp_solver = AcadosOcpSolver(ocp);

x_traj_init = [linspace(0, 0, N+1); linspace(pi, 0, N+1); linspace(0, 0, N+1); linspace(0, 0, N+1)];
u_traj_init = zeros(nu, N);

% if not set, the trajectory is initialized with the previous solution
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);

% solve
tic;

matlab_sqp_loop = 0;
if ~matlab_sqp_loop
    % solve ocp
    ocp_solver.solve();
else

    % do one step at a time
    ocp_solver.set('nlp_solver_max_iter', 1);

    for ii=1:nlp_solver_max_iter

        disp(['iteration number ', num2str(ii)])

        % solve the system using 1 SQP iteration
        ocp_solver.solve();

        % print 1-iteration stat
        ocp_solver.print('stat');

        % check stability of qp
        qp_A = ocp_solver.get('qp_A');
        qp_A_eig_max = 0;
        for jj=1:length(qp_A)
            tmp_A = qp_A{jj};
            qp_A_eig = eig(tmp_A);
            tmp = max(abs(qp_A_eig));
            if tmp>qp_A_eig_max
                qp_A_eig_max = tmp;
            end
        end
        fprintf('A eig max %e\n', qp_A_eig_max);

        % compute condition number and eigenvalues of hessian of (partial) cond qp
        qp_cond_H = ocp_solver.get('qp_solver_cond_H');
        if iscell(qp_cond_H)

            for jj=1:length(qp_cond_H)

                tmp_H = qp_cond_H{jj};
                nv = size(tmp_H, 1);
                % make full
                for jj=1:nv
                    for ii=jj+1:nv
                        tmp_H(jj,ii) = tmp_H(ii,jj);
                    end
                end
                qp_H_cond_num = cond(tmp_H);
                qp_H_eig = eig(tmp_H);
                fprintf('cond H condition number %e min eigenval %e max eigenvalue %e\n', qp_H_cond_num, min(qp_H_eig), max(qp_H_eig));

            end
        else
            nv = size(qp_cond_H, 1);
            % make full
            for jj=1:nv
                for ii=jj+1:nv
                    qp_cond_H(jj,ii) = qp_cond_H(ii,jj);
                end
            end
            qp_H_cond_num = cond(qp_cond_H);
            qp_H_eig = eig(qp_cond_H);
            fprintf('cond H condition number %e min eigenval %e max eigenvalue %e\n', qp_H_cond_num, min(qp_H_eig), max(qp_H_eig));

        end

		% check residuals and terminate if tol is reached
		residuals = ocp_solver.get('residuals');
		if residuals(1) < nlp_solver_tol_stat && residuals(2) < nlp_solver_tol_eq && residuals(3) < nlp_solver_tol_ineq && residuals(4) < nlp_solver_tol_comp
			break
		end
    end
end

time_ext = toc;

% get solution
u = ocp_solver.get('u');
x = ocp_solver.get('x');

%% evaluation
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_sim = ocp_solver.get('time_sim');
time_lin_remaining = time_lin - time_sim; % linearization excluding integrators
time_reg = ocp_solver.get('time_reg');
time_qp_solver_call = ocp_solver.get('time_qp_solver_call');
time_qp_sol = ocp_solver.get('time_qp_sol');
time_qp_remaining = time_qp_sol - time_qp_solver_call; % time for QP solution in addition to QP solver call, mostly condensing
time_remaining_nlp = time_tot - time_lin - time_reg - time_qp_sol;

fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

ocp_solver.print('stat');

%%
figure;
bar_data = 1e3*[time_sim, time_lin_remaining, time_reg, time_qp_solver_call, time_qp_remaining, time_remaining_nlp]';
bar(0, bar_data, 'stacked')
legend('integrators', 'remaining linearization', 'regularization', 'QP solver call', 'QP processing: condensing, etc', 'remaining')
xlim([-.5, 1])
ylabel('time in [ms]');
xticks([]);

%% plot trajectory
figure;
subplot(2,1,1);
plot(0:N, x);
xlim([0 N]);
legend('p', 'theta', 'v', 'omega');
subplot(2,1,2);
plot(0:N-1, u);
xlim([0 N]);
legend('F');


%% plot residual
% stat = ocp_solver.get('stat');
% if (strcmp(nlp_solver, 'sqp'))
%     figure;
%     plot([0: size(stat,1)-1], log10(stat(:,2)), 'r-x');
%     hold on
%     plot([0: size(stat,1)-1], log10(stat(:,3)), 'b-x');
%     plot([0: size(stat,1)-1], log10(stat(:,4)), 'g-x');
%     plot([0: size(stat,1)-1], log10(stat(:,5)), 'k-x');
% %    semilogy(0: size(stat,1)-1, stat(:,2), 'r-x');
% %    hold on
% %    semilogy(0: size(stat,1)-1, stat(:,3), 'b-x');
% %    semilogy(0: size(stat,1)-1, stat(:,4), 'g-x');
% %    semilogy(0: size(stat,1)-1, stat(:,5), 'k-x');
%     hold off
%     xlabel('iter')
%     ylabel('res')
%     legend('res stat', 'res eq', 'res ineq', 'res compl');
% end


if status==0
    fprintf('\nsuccess!\n\n');
else
    fprintf('\nsolution failed!\n\n');
end


if is_octave()
    waitforbuttonpress;
end

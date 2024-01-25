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

import casadi.*
addpath('../pendulum_on_cart_model/');

for itest = 1:3
    %% arguments
    compile_interface = 'auto';
    codgen_model = 'true';
    gnsf_detect_struct = 'true';

    % discretization
    N = 100;
    h = 0.01;

    nlp_solver = 'sqp';
    %nlp_solver = 'sqp_rti';
    nlp_solver_exact_hessian = 'false';
    %nlp_solver_exact_hessian = 'true';
    regularize_method = 'no_regularize';
    %regularize_method = 'project';
    %regularize_method = 'project_reduc_hess';
    %regularize_method = 'mirror';
    %regularize_method = 'convexify';
    nlp_solver_max_iter = 100;
    test_tol = 2e-8;
    nlp_solver_tol = 1e-10;
    nlp_solver_tol_stat = nlp_solver_tol;
    nlp_solver_tol_eq   = nlp_solver_tol;
    nlp_solver_tol_ineq = nlp_solver_tol;
    nlp_solver_tol_comp = nlp_solver_tol;
    nlp_solver_ext_qp_res = 1;
    %qp_solver = 'partial_condensing_hpipm';
    %qp_solver = 'full_condensing_hpipm';
    qp_solver = 'full_condensing_qpoases';
    qp_solver_cond_N = 5;
    qp_solver_cond_ric_alg = 0;
    qp_solver_ric_alg = 0;
    qp_solver_warm_start = 1;
    qp_solver_iter_max = 100;
    %sim_method = 'erk';
    sim_method = 'irk';
%     sim_method = 'irk_gnsf';
    sim_method_num_stages = 4;
    sim_method_num_steps = 3;

    if itest == 1
        cost_type = 'linear_ls';
    elseif itest == 2
        cost_type = 'ext_cost';
    else
        cost_type = 'auto';
    end
    model_name = ['pendulum_' num2str(itest)];

    %% create model entries
    model = pendulum_on_cart_model();

    % dims
    T = N*h; % horizon length time
    nx = model.nx;
    nu = model.nu;
    ny = nu+nx; % number of outputs in lagrange term
    ny_e = nx; % number of outputs in mayer term

    nbx = 0;
    nbu = 0;

    % cost
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

    % constraints
    x0 = [0; pi; 0; 0];
    %Jbx = zeros(nbx, nx); for ii=1:nbx Jbx(ii,ii)=1.0; end
    %lbx = -4*ones(nbx, 1);
    %ubx =  4*ones(nbx, 1);
    Jbu = zeros(nbu, nu); for ii=1:nbu Jbu(ii,ii)=1.0; end
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
    if (strcmp(cost_type, 'linear_ls'))
        ocp_model.set('cost_Vu', Vu);
        ocp_model.set('cost_Vx', Vx);
        ocp_model.set('cost_Vx_e', Vx_e);
        ocp_model.set('cost_W', W);
        ocp_model.set('cost_W_e', W_e);
        ocp_model.set('cost_y_ref', yr);
        ocp_model.set('cost_y_ref_e', yr_e);
    else % if (strcmp(cost_type, 'ext_cost'))
        ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
        ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
    end
    % dynamics
    if (strcmp(sim_method, 'erk'))
        ocp_model.set('dyn_type', 'explicit');
        ocp_model.set('dyn_expr_f', model.expr_f_expl);
    else % irk irk_gnsf
        ocp_model.set('dyn_type', 'implicit');
        ocp_model.set('dyn_expr_f', model.expr_f_impl);
    end

    %% constraints
    ocp_model.set('constr_x0', x0);
    if itest == 1
        ng = 0;
        ocp_model.set('constr_expr_h', model.expr_h);
        ocp_model.set('constr_lh', lbu);
        ocp_model.set('constr_uh', ubu);
        ocp_model.set('constr_expr_h_0', model.expr_h);
        ocp_model.set('constr_lh_0', lbu);
        ocp_model.set('constr_uh_0', ubu);
    elseif itest == 2
        ng = 1;
        C = zeros(ng, nx);
        D = zeros(ng, nu);
        D(1, nu) = 1;
        ocp_model.set('constr_D', D);
        ocp_model.set('constr_C', C);

        ocp_model.set('constr_lg', lbu);
        ocp_model.set('constr_ug', ubu);

    elseif itest == 3
        ng = 0;
        ocp_model.set('constr_type', 'auto');
        ocp_model.set('constr_expr_h', model.expr_h);
        ocp_model.set('constr_lh', lbu);
        ocp_model.set('constr_uh', ubu);
%         ocp_model.set('constr_type_e', 'auto');
%         ocp_model.set('constr_expr_h_e', SX.sym('terminal_constraint',0,0));
%         ocp_model.set('constr_lh_e', []);
%         ocp_model.set('constr_uh_e', []);
    end

    %% acados ocp opts
    ocp_opts = acados_ocp_opts();
    ocp_opts.set('compile_interface', compile_interface);
    ocp_opts.set('codgen_model', codgen_model);
    ocp_opts.set('param_scheme_N', N);
    ocp_opts.set('nlp_solver', nlp_solver);
    ocp_opts.set('nlp_solver_exact_hessian', nlp_solver_exact_hessian);
    % ocp_opts.set('nlp_solver_exact_hessian', 'TEST');
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
    ocp_opts.set('qp_solver_iter_max', qp_solver_iter_max);
    ocp_opts.set('sim_method', sim_method);
    ocp_opts.set('sim_method_num_stages', sim_method_num_stages);
    ocp_opts.set('sim_method_num_steps', sim_method_num_steps);
    if (strcmp(sim_method, 'irk_gnsf'))
        ocp_opts.set('gnsf_detect_struct', gnsf_detect_struct);
    end
    ocp_opts.set('ext_fun_compile_flags', '');

    %% acados ocp
    % create ocp
    ocp = acados_ocp(ocp_model, ocp_opts);

    % set trajectory initialization
    %x_traj_init = zeros(nx, N+1);
    %for ii=1:N x_traj_init(:,ii) = [0; pi; 0; 0]; end
    x_traj_init = [linspace(0, 0, N+1); linspace(pi, 0, N+1); linspace(0, 0, N+1); linspace(0, 0, N+1)];

    u_traj_init = zeros(nu, N);
    ocp.set('init_x', x_traj_init);
    ocp.set('init_u', u_traj_init);

    % change number of sqp iterations
    ocp.set('nlp_solver_max_iter', 20);

    % modify numerical data for a certain stage
    some_stages = 1:10:N-1;
    for i = some_stages
        if (strcmp(cost_type, 'linear_ls'))
            ocp.set('cost_Vx', Vx, i); % cost_y_ref, cost_Vu, cost_Vx, cost_W, cost_Z, cost_Zl,...
             % cost_Zu, cost_z, cost_zl, cost_zu;
            ocp.set('cost_Vu', Vu, i);
            ocp.set('cost_y_ref', yr, i);
        end
        if ng > 0
            ocp.set('constr_C', C, i);
            ocp.set('constr_D', D, i);
            ocp.set('constr_ug', ubu, i);
            ocp.set('constr_lg', lbu, i);
        end
    end

    % solve
    tic;
    ocp.solve();
    time_ext=toc;
    % get solution
    utraj = ocp.get('u');
    xtraj = ocp.get('x');

    %% evaluation
    status = ocp.get('status');
    sqp_iter = ocp.get('sqp_iter');
    time_tot = ocp.get('time_tot');
    time_lin = ocp.get('time_lin');
    time_reg = ocp.get('time_reg');
    time_qp_sol = ocp.get('time_qp_sol');

    fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

    stat = ocp.get('stat');
    ocp.print('stat')

    if itest == 1
        utraj_ref = utraj;
        xtraj_ref = xtraj;
    else
        err_x = max(max(abs(xtraj - xtraj_ref)));
        err_u = max(max(abs(utraj - utraj_ref)));
        if max(err_x, err_u) > test_tol
            error(['\nSolutions differ by more than test_tol = ' num2str(test_tol) ' : ' num2str(err_x) ' , ' num2str(err_u)]);
        end
    end

    if status~=0
        error('test_ocp_pendulum_on_cart: solution failed!');
    elseif test_tol < max(stat(end,2:5))
        error('test_ocp_pendulum_on_cart: residuals bigger than test_tol!');
    elseif sqp_iter > 11
        error('test_ocp_pendulum_on_cart: sqp_iter > 11, this problem is typically solved within less iterations!');
    end
% For debugging
%     figure;
%     plot(1:N+1, xtraj);
%     legend('p', 'theta', 'v', 'omega');

end
fprintf('\ntest_ocp_pendulum_on_cart: success!\n');

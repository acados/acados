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


import casadi.*
addpath('../pendulum_on_cart_model/');
check_acados_requirements();

for itest = 1:3
    %% arguments
    N = 100;
    h = 0.01;
    T = N*h;

    test_tol = 2e-8;
    sim_method = 'IRK';
    sim_method_num_stages = 4;
    sim_method_num_steps = 3;

    if itest == 1
        cost_type = 'LINEAR_LS';
    elseif itest == 2
        cost_type = 'EXTERNAL';
    else
        cost_type = 'AUTO';
    end
    model_name = ['pendulum_' num2str(itest)];

    %% create model entries
    model = pendulum_on_cart_model();
    nx = model.nx;
    nu = model.nu;
    ny = nu + nx;
    ny_e = nx;

    Vu = zeros(ny, nu);
    for ii = 1:nu
        Vu(ii, ii) = 1.0;
    end
    Vx = zeros(ny, nx);
    for ii = 1:nx
        Vx(nu + ii, ii) = 1.0;
    end
    Vx_e = zeros(ny_e, nx);
    for ii = 1:nx
        Vx_e(ii, ii) = 1.0;
    end

    W = eye(ny);
    for ii = 1:nu
        W(ii, ii) = 1e-2;
    end
    for ii = nu + 1:nu + nx/2
        W(ii, ii) = 1e3;
    end
    for ii = nu + nx/2 + 1:nu + nx
        W(ii, ii) = 1e-2;
    end
    W_e = W(nu + 1:nu + nx, nu + 1:nu + nx);
    yr = zeros(ny, 1);
    yr_e = zeros(ny_e, 1);

    x0 = [0; pi; 0; 0];
    lbu = -80 * ones(nu, 1);
    ubu =  80 * ones(nu, 1);

    %% acados OCP model
    acados_model = AcadosModel();
    acados_model.name = model_name;
    acados_model.x = model.sym_x;
    if isfield(model, 'sym_u')
        acados_model.u = model.sym_u;
    end
    if isfield(model, 'sym_xdot')
        acados_model.xdot = model.sym_xdot;
    end

    if strcmp(sim_method, 'ERK')
        acados_model.f_expl_expr = model.dyn_expr_f_expl;
    else
        acados_model.f_impl_expr = model.dyn_expr_f_impl;
    end

    ocp = AcadosOcp();
    ocp.model = acados_model;
    ocp.solver_options.N_horizon = N;
    ocp.solver_options.tf = T;
    ocp.solver_options.nlp_solver_type = 'SQP';
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    ocp.solver_options.integrator_type = sim_method;
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES';
    ocp.solver_options.qp_solver_cond_N = 5;
    ocp.solver_options.qp_solver_warm_start = 1;
    ocp.solver_options.qp_solver_iter_max = 100;
    ocp.solver_options.nlp_solver_max_iter = 30;
    ocp.solver_options.nlp_solver_tol_stat = 1e-10;
    ocp.solver_options.nlp_solver_tol_eq = 1e-10;
    ocp.solver_options.nlp_solver_tol_ineq = 1e-10;
    ocp.solver_options.nlp_solver_tol_comp = 1e-10;
    ocp.solver_options.nlp_solver_ext_qp_res = 1;
    ocp.solver_options.sim_method_num_stages = sim_method_num_stages;
    ocp.solver_options.sim_method_num_steps = sim_method_num_steps;
    ocp.solver_options.regularize_method = 'NO_REGULARIZE';

    %% cost
    ocp.cost.cost_type = cost_type;
    ocp.cost.cost_type_e = cost_type;
    if strcmp(cost_type, 'LINEAR_LS')
        ocp.cost.Vu = Vu;
        ocp.cost.Vx = Vx;
        ocp.cost.Vx_e = Vx_e;
        ocp.cost.W = W;
        ocp.cost.W_e = W_e;
        ocp.cost.yref = yr;
        ocp.cost.yref_e = yr_e;
    else
        ocp.model.cost_expr_ext_cost = model.cost_expr_ext_cost;
        ocp.model.cost_expr_ext_cost_e = model.cost_expr_ext_cost_e;
    end

    %% constraints
    ocp.constraints.x0 = x0;
    if itest == 1
        ocp.model.con_h_expr = model.constr_expr_h;
        ocp.model.con_h_expr_0 = model.constr_expr_h;
        ocp.constraints.lh = lbu;
        ocp.constraints.uh = ubu;
        ocp.constraints.lh_0 = lbu;
        ocp.constraints.uh_0 = ubu;
    elseif itest == 2
        ocp.model.con_h_expr = model.constr_expr_h;
        ocp.model.con_h_expr_0 = model.constr_expr_h;
        ocp.constraints.lh = lbu;
        ocp.constraints.uh = ubu;
        ocp.constraints.lh_0 = lbu;
        ocp.constraints.uh_0 = ubu;
    else
        ocp.model.con_h_expr = model.constr_expr_h;
        ocp.model.con_h_expr_0 = model.constr_expr_h;
        ocp.constraints.lh = lbu;
        ocp.constraints.uh = ubu;
        ocp.constraints.lh_0 = lbu;
        ocp.constraints.uh_0 = ubu;
    end

    %% acados OCP solver
    ocp_solver = AcadosOcpSolver(ocp);

    % set trajectory initialization
    x_traj_init = [linspace(0, 0, N + 1); linspace(pi, 0, N + 1); linspace(0, 0, N + 1); linspace(0, 0, N + 1)];
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

    %% evaluation
    status = ocp_solver.get('status');
    sqp_iter = ocp_solver.get('sqp_iter');
    time_tot = ocp_solver.get('time_tot');
    time_lin = ocp_solver.get('time_lin');
    time_reg = ocp_solver.get('time_reg');
    time_qp_sol = ocp_solver.get('time_qp_sol');

    fprintf('\nstatus = %d, sqp_iter = %d, time_ext = %f [ms], time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n', ...
        status, sqp_iter, time_ext*1e3, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3);

    stat = ocp_solver.get('stat');
    ocp_solver.print('stat');

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

    if status ~= 0
        error('test_ocp_pendulum_on_cart: solution failed!');
    elseif test_tol < max(stat(end, 2:5))
        error('test_ocp_pendulum_on_cart: residuals bigger than test_tol!');
    elseif sqp_iter > 11
        error('test_ocp_pendulum_on_cart: sqp_iter > 11, this problem is typically solved within less iterations!');
    end

    % For debugging
    % figure;
    % plot(1:N + 1, xtraj);
    % legend('p', 'theta', 'v', 'omega');
end

fprintf('\ntest_ocp_pendulum_on_cart: success!\n');

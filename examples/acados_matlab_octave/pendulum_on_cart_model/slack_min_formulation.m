% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


function xtraj = slack_min_formulation(formulation)
    if nargin < 1
        formulation = 's_slack';
    end

    import casadi.*

    ACADOS_INFTY = get_acados_infty();

    % create ocp object to formulate the OCP
    ocp = AcadosOcp();

    % set model
    model = get_pendulum_on_cart_model();
    ocp.model = model;
    ocp.model.name = formulation;

    Tf = 1.0;
    N = 10;

    % set prediction horizon
    ocp.solver_options.N_horizon = N;
    ocp.solver_options.tf = Tf;

    % cost matrices
    Q_mat = 2*diag([1e3, 1e3, 1e-2, 1e-2]);
    R_mat = 2*diag([1e-2]);

    % path cost
    ocp.cost.cost_type = 'EXTERNAL';
    W_mat = blkdiag(Q_mat, R_mat);
    ocp.model.cost_expr_ext_cost = 0.5 * (vertcat(model.x, model.u))' * W_mat * vertcat(model.x, model.u);
    % terminal cost
    ocp.cost.cost_type_e = 'EXTERNAL';
    ocp.model.cost_expr_ext_cost_e = 0.5 * model.x' * Q_mat * model.x;

    % set constraints
    Fmax = 40;
    ocp.constraints.lbu = [-Fmax];
    ocp.constraints.ubu = [+Fmax];
    ocp.constraints.idxbu = 0;

    % initial condition
    ocp.constraints.x0 = [0.0; 0.2*pi; 0.0; 0.0];

    fprintf('using formulation %s\n', formulation);

    % add cost term min(x[0], x[2]) to cost via slack: s <= x[0], s <= x[2]
    if strcmp(formulation, 'u_slack')
        % add u
        new_u = SX.sym('new_u', 1, 1);
        ocp.model.u = vertcat(model.u, new_u);
        % add constraints u <= x_i <=> -inf <= u - x_i <= 0
        ocp.model.con_h_expr = vertcat(new_u - model.x(1), new_u - model.x(4));
        ocp.constraints.uh = zeros(2, 1);
        ocp.constraints.lh = -ACADOS_INFTY * ones(2, 1);

        ocp.model.con_h_expr_0 = ocp.model.con_h_expr;
        ocp.constraints.uh_0 = ocp.constraints.uh;
        ocp.constraints.lh_0 = ocp.constraints.lh;

        % add cost -u
        ocp.model.cost_expr_ext_cost = ocp.model.cost_expr_ext_cost - new_u;
    elseif strcmp(formulation, 'u_slack2')
        % add u
        new_u = SX.sym('new_u', 1, 1);
        ocp.model.u = vertcat(model.u, new_u);
        % add constraints u <= x_i <=> -inf <= u - x_i <= 0
        ocp.model.con_h_expr = vertcat(new_u + model.x(1), new_u + model.x(4));
        ocp.constraints.uh = ACADOS_INFTY * ones(2, 1);
        ocp.constraints.lh = zeros(2, 1);

        ocp.model.con_h_expr_0 = ocp.model.con_h_expr;
        ocp.constraints.uh_0 = ocp.constraints.uh;
        ocp.constraints.lh_0 = ocp.constraints.lh;

        % add cost u
        ocp.model.cost_expr_ext_cost = ocp.model.cost_expr_ext_cost + new_u;
    elseif strcmp(formulation, 's_slack')
        % add s
        ns = 1;
        % add constraints: s <= x_i
        ocp.model.con_h_expr = vertcat(model.x(1), model.x(4));
        ocp.constraints.uh = ACADOS_INFTY * ones(2, 1);
        ocp.constraints.lh = zeros(2, 1);
        ocp.constraints.idxs_rev = [-1; 0; 0];
        ocp.constraints.ls = -ACADOS_INFTY * ones(ns, 1);
        ocp.constraints.us = zeros(ns, 1);
        ocp.cost.zl = 1.0;
        ocp.cost.Zl = -0.0;
        ocp.cost.zu = 1.0;
        ocp.cost.Zu = 0.0;

        ocp.model.con_h_expr_0 = ocp.model.con_h_expr;
        ocp.constraints.uh_0 = ocp.constraints.uh;
        ocp.constraints.lh_0 = ocp.constraints.lh;
        ocp.cost.zl_0 = ocp.cost.zl;
        ocp.cost.Zl_0 = ocp.cost.Zl;
        ocp.cost.zu_0 = ocp.cost.zu;
        ocp.cost.Zu_0 = ocp.cost.Zu;
        nbx_0 = numel(ocp.constraints.lbx_0);
        nbu = numel(ocp.constraints.lbu);
        ocp.constraints.idxs_rev_0 = [-ones(nbx_0+nbu,1); 0; 0];
        ocp.constraints.ls_0 = ocp.constraints.ls;
        ocp.constraints.us_0 = ocp.constraints.us;
    end
    ocp.solver_options.qp_solver_t0_init = 0;

    % set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    ocp.solver_options.integrator_type = 'IRK';
    ocp.solver_options.nlp_solver_type = 'SQP';
    % ocp.solver_options.print_level = 5;
    % ocp.solver_options.nlp_solver_max_iter = 2;

    nx = length(model.x);
    nu = length(model.u);

    ocp_solver = AcadosOcpSolver(ocp);

    xtraj = zeros(N+1, nx);
    utraj = zeros(N, nu);

    ocp_solver.solve();
    status = ocp_solver.get('status');
    ocp_solver.print_statistics();

    if status ~= 0
        error('acados returned status %d.', status);
    end

    % get solution
    for i = 1:N
        xtraj(i,:) = ocp_solver.get('x', i-1);
        utraj(i,:) = ocp_solver.get('u', i-1);
    end
    xtraj(N+1,:) = ocp_solver.get('x', N);

    min_x_vals = min(xtraj(:,1), xtraj(:,4));
    if strcmp(formulation, 'u_slack')
        slack_vals = utraj(:,2);
        assert(all(abs(min_x_vals(1:end-1) - slack_vals) < 1e-6));
    elseif strcmp(formulation, 'u_slack2')
        slack_vals = utraj(:,2);
        assert(all(abs(min_x_vals(1:end-1) + slack_vals) < 1e-6));
    elseif strcmp(formulation, 's_slack')
        slack_vals = zeros(N, 1);
        unused_slack_vals = zeros(N, 1);
        for i = 1:N
            sl = ocp_solver.get('sl', i-1);
            su = ocp_solver.get('su', i-1);
            slack_vals(i) = sl(1);
            unused_slack_vals(i) = su(1);
        end
        assert(all(abs(min_x_vals(1:end-1) + slack_vals) < 1e-6));
        disp('unused_slack_vals=');
        disp(unused_slack_vals);
        % plot slacks
        utraj = [utraj, slack_vals];
    end
end


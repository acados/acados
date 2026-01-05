function solver = setup_ocp_solver(x0, umax, dt_0, N_horizon, Tf, tol, with_abs_cost, hessian_approx, regularize_method, anderson_activation_threshold)

    model = get_furuta_model();

    nx = length(model.x);
    nu = 1;
    ny = nx + nu;
    ny_e = nx;

    ocp = AcadosOcp();
    ocp.model = model;
    ocp.solver_options.N_horizon = N_horizon;

    % cost
    ocp.cost.cost_type = 'NONLINEAR_LS';
    ocp.cost.cost_type_e = 'NONLINEAR_LS';

    Q_mat = diag([50.0, 500.0, 1.0, 1.0]);
    R_mat = diag([1e3]);

    ocp.cost.W   = blkdiag(Q_mat, R_mat);
    ocp.cost.W_e = Q_mat;

    ocp.model.cost_y_expr   = [model.x; model.u];
    ocp.model.cost_y_expr_e = model.x;
    ocp.cost.yref   = zeros(ny,1);
    ocp.cost.yref_e = zeros(ny_e,1);

    % constraints
    ocp.constraints.lbu   = -umax;
    ocp.constraints.ubu   = +umax;
    ocp.constraints.idxbu = 0;

    if with_abs_cost
        val = 1.4;
        % add terminal / stage box constraint and slacks similar to python example
        ocp.constraints.idxbx = 0;
        ocp.constraints.lbx = val;
        ocp.constraints.ubx = val;
        ocp.constraints.idxsbx = 0;

        ocp.cost.zl = 1e3 * 1.0;
        ocp.cost.zu = 1e3 * 1.0;
        ocp.cost.Zl = 0.0 * 1.0;
        ocp.cost.Zu = 0.0 * 1.0;
    end

    ocp.constraints.x0 = x0;

    % solver options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = hessian_approx;
    ocp.solver_options.regularize_method = regularize_method;
    ocp.solver_options.integrator_type = 'ERK';
    ocp.solver_options.reg_epsilon = 5e-2;

    % non-uniform grid
    ocp.solver_options.time_steps = [dt_0, repmat((Tf-dt_0)/(N_horizon-1), 1, N_horizon-1)];
    ocp.solver_options.sim_method_num_steps = [1, repmat(2,1,N_horizon-1)];
    ocp.solver_options.levenberg_marquardt = 1e-6;
    ocp.solver_options.with_anderson_acceleration = true;
    ocp.solver_options.anderson_activation_threshold = anderson_activation_threshold;

    ocp.solver_options.qp_solver_cond_N = N_horizon;
    ocp.solver_options.nlp_solver_tol_stat = tol;
    ocp.solver_options.nlp_solver_tol_eq = tol;
    ocp.solver_options.nlp_solver_tol_ineq = tol;
    ocp.solver_options.nlp_solver_tol_comp = tol;
    ocp.solver_options.tf = Tf;

    solver = AcadosOcpSolver(ocp);

end

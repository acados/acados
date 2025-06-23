function test_qpscaling()
    import casadi.*

    % solver name
    solver_name = 'test_qpscaling';

    % create solver
    nlp_solver_type = 'SQP';
    use_qp_scaling = true;
    soft_h = true;
    ocp_solver = create_solver(solver_name, nlp_solver_type, use_qp_scaling, soft_h);

    % solve
    ocp_solver.solve();

    if use_qp_scaling
        constraints_scaling = ocp_solver.get_qp_scaling_constraints(0);
        objective_scaling = ocp_solver.get_qp_scaling_objective();
        qpscaling_status = ocp_solver.get('qpscaling_status');
        fprintf('constraints scaling: %f\n', constraints_scaling);
        fprintf('objective scaling: %f\n', objective_scaling);
        fprintf('QP scaling status: %d\n', qpscaling_status);
    end

    % print solution
    ocp_solver.print('stat');
end


function ocp_solver = create_solver(solver_name, nlp_solver_type, use_qp_scaling, soft_h)
    import casadi.*

    ACADOS_INFTY = get_acados_infty();
    ocp = AcadosOcp();

    nx = 2;
    % set model
    ocp.model.name = ['dense_nlp_' solver_name];
    ocp.model.x = SX.sym('x', nx, 1);

    ny = nx;

    % discretization
    N = 0;

    ocp.cost.W_e = 1e3*2*diag([1e3, 1e3]);

    ocp.cost.cost_type = 'LINEAR_LS';
    ocp.cost.cost_type_e = 'LINEAR_LS';

    ocp.cost.Vx_e = eye(nx);
    ocp.cost.yref_e = ones(ny, 1);

    % set constraints
    xmax = 2.0;
    ocp.constraints.lbx_e = -xmax * ones(nx,1);
    ocp.constraints.ubx_e = xmax * ones(nx,1);
    ocp.constraints.idxbx_e = 0:(nx-1);

    % define soft nonlinear constraint
    scale_h = 1.0;
    radius = 1.0;
    ocp.model.con_h_expr_e = scale_h * (ocp.model.x(1)^2 + ocp.model.x(2)^2);
    ocp.constraints.lh_e = -1000 * ones(1,1);
    ocp.constraints.lh_e = -ACADOS_INFTY * ones(1,1);
    ocp.constraints.uh_e = scale_h * radius^2 * ones(1,1);

    % soften
    if soft_h
        ocp.constraints.idxsh_e = 0;
        ocp.cost.zl_e = [ocp.cost.zl_e; 1.0];
        ocp.cost.zu_e = [ocp.cost.zu_e; 1.0];
        ocp.cost.Zl_e = [ocp.cost.Zl_e; 1.0];
        ocp.cost.Zu_e = [ocp.cost.Zu_e; 1.0];
    end

    % set options
    solver_options = ocp.solver_options;
    solver_options.N_horizon = N;

    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    qp_tol = 5e-9;
    solver_options.qp_solver_tol_stat = qp_tol;
    solver_options.qp_solver_tol_eq = qp_tol;
    solver_options.qp_solver_tol_ineq = qp_tol;
    solver_options.qp_solver_tol_comp = qp_tol;
    solver_options.qp_solver_ric_alg = 1;
    solver_options.qp_solver_mu0 = 1e4;
    solver_options.qp_solver_iter_max = 400;
    solver_options.hessian_approx = 'GAUSS_NEWTON';
    solver_options.nlp_solver_type = nlp_solver_type;
    solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH';
    solver_options.globalization_full_step_dual = true;
    % solver_options.print_level = 1;
    % solver_options.nlp_solver_max_iter = 2;
    solver_options.nlp_solver_ext_qp_res = 0;

    if use_qp_scaling
        ocp.solver_options.qpscaling_scale_constraints = 'INF_NORM';
        ocp.solver_options.qpscaling_scale_objective = 'OBJECTIVE_GERSHGORIN';
    end

    % create ocp solver
    ocp_solver = AcadosOcpSolver(ocp);
end




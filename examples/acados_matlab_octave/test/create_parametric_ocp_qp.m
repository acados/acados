function [ocp, x0] = create_parametric_ocp_qp(N, np)
    import casadi.*

    nx = 3; % state size
    nu = 3; % input size
    ny = nx+nu;
    ny_e = nx;

    x0 = ones(nx, 1); % initial state

    %% model
    model = AcadosModel();

    sym_x = SX.sym('x', nx);
    sym_u = SX.sym('u', nu);
    sym_p = SX.sym('p', np);

    model.x = sym_x;
    model.u = sym_u;
    model.p = sym_p;
    model.disc_dyn_expr = sym_x + sym_p(1);

    model.name = 'trivial';

    ocp = AcadosOcp();
    ocp.model = model;

    % cost
    cost_type = 'LINEAR_LS';
    cost_type_e = 'LINEAR_LS';
    Vx = zeros(ny,nx); Vx(1:nx,:) = eye(nx);        % state-to-output matrix in lagrange term
    Vu = zeros(ny,nu); Vu(nx+1:ny,:) = eye(nu);     % input-to-output matrix in lagrange term
    Vx_e = zeros(ny_e,nx); Vx_e(1:nx,:) = eye(nx);  % state-to-output matrix in mayer term
    W = eye(ny);
    W_e = 5 * W(1:ny_e,1:ny_e);                         % cost weights in mayer term
    y_ref = zeros(ny,1);                            % set initial references
    y_ref_e = zeros(ny_e,1);

    ocp.cost.cost_type = cost_type;
    ocp.cost.cost_type_e = cost_type_e;
    ocp.cost.Vx = Vx;
    ocp.cost.Vu = Vu;
    ocp.cost.Vx_e = Vx_e;
    ocp.cost.W = W;
    ocp.cost.W_e = W_e;
    ocp.cost.yref = y_ref;
    ocp.cost.yref_e = y_ref_e;

    % constraints
    ocp.constraints.x0 = x0;  % set the initial state

    %% options
    ocp.solver_options.N_horizon = N;
    ocp.solver_options.tf = 1;
    ocp.solver_options.nlp_solver_type = 'SQP';
    ocp.solver_options.integrator_type = 'DISCRETE';

    %% Simulink opts
    simulink_opts = get_acados_simulink_opts;

    % inputs
    simulink_opts.inputs.y_ref = 0;
    simulink_opts.inputs.y_ref_0 = 0;
    simulink_opts.inputs.y_ref_e = 0;

    simulink_opts.inputs.cost_W = 0;
    simulink_opts.inputs.cost_W_0 = 0;
    simulink_opts.inputs.cost_W_e = 0;

    simulink_opts.inputs.x_init = 1;
    simulink_opts.inputs.u_init = 1;
    simulink_opts.inputs.pi_init = 1;

    simulink_opts.inputs.reset_solver = 1;
    simulink_opts.inputs.ignore_inits = 1;

    % outputs
    simulink_opts.outputs.u0 = 0;
    simulink_opts.outputs.utraj = 1;
    simulink_opts.outputs.xtraj = 1;
    simulink_opts.outputs.pi_all = 1;

    simulink_opts.outputs.solver_status = 1;
    simulink_opts.outputs.CPU_time = 0;
    simulink_opts.outputs.sqp_iter = 1;
    simulink_opts.outputs.x1 = 0;

    ocp.simulink_opts = simulink_opts;
end

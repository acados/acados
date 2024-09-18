
function ocp = create_ocp_formulation(p_global, m, l, C, lut, use_p_global, p_global_values)
    ocp = AcadosOcp();

    % Set model
    model = export_pendulum_ode_model(p_global, m, l, C, lut);
    model.p_global = p_global;
    ocp.model = model;

    % Dimensions
    nx = size(model.x, 1);
    nu = size(model.u, 1);
    ny = nx + nu;
    ny_e = nx;

    % Cost settings
    Q = diag([1e3, 1e3, 1e-2, 1e-2]);
    R = diag([1e-2]);

    ocp.cost.W_e = Q;
    ocp.cost.W = blkdiag(Q, R);

    ocp.cost.cost_type = 'LINEAR_LS';
    ocp.cost.cost_type_e = 'LINEAR_LS';

    ocp.cost.Vx = zeros(ny, nx);
    ocp.cost.Vx(1:nx, 1:nx) = eye(nx);

    Vu = zeros(ny, nu);
    Vu(5, 1) = 1.0;
    ocp.cost.Vu = Vu;

    ocp.cost.Vx_e = eye(nx);

    ocp.cost.yref = zeros(ny, 1);
    ocp.cost.yref_e = zeros(ny_e, 1);

    % Constraints
    Fmax = 80;
    ocp.constraints.lbu = -Fmax;
    ocp.constraints.ubu = Fmax;
    ocp.constraints.idxbu = 0;

    ocp.constraints.x0 = [0.0; pi; 0.0; 0.0];

    % Set options
    Tf = 1.0;
    N_horizon = 20;

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    ocp.solver_options.integrator_type = 'ERK';
    ocp.solver_options.print_level = 0;
    ocp.solver_options.nlp_solver_type = 'SQP_RTI';

    ocp.solver_options.tf = Tf;
    ocp.solver_options.N_horizon = N_horizon;

    % Parameters
    ocp.parameter_values = 9.81;

    if ~use_p_global
        model.p = vertcat(model.p, p_global);
        model.p_global = [];
        ocp.parameter_values = [ocp.parameter_values; p_global_values];
    end
end


function model = export_pendulum_ode_model(p_global, m, l, C, lut)
    import casadi.*
    model_name = 'pendulum';

    % Constants
    m_cart = 1.0; % mass of the cart [kg]

    % Parameters
    g = MX.sym('g');
    p = g;

    % States & controls
    x1 = MX.sym('x1');
    theta = MX.sym('theta');
    v1 = MX.sym('v1');
    dtheta = MX.sym('dtheta');
    x = vertcat(x1, theta, v1, dtheta);

    F = MX.sym('F');
    u = F;

    xdot = MX.sym('xdot', length(x));

    % Dynamics
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    denominator = m_cart + m - m*cos_theta^2;
    f_expl = vertcat(v1, dtheta, ...
                     (-m*l*sin_theta*dtheta^2 + m*g*cos_theta*sin_theta + F) / denominator, ...
                     (-m*l*cos_theta*sin_theta*dtheta^2 + F*cos_theta + (m_cart + m)*g*sin_theta) / (l*denominator));

    if lut
        knots = {[0,0,0,0,0.2,0.5,0.8,1,1,1,1],[0,0,0,0.1,0.5,0.9,1,1,1]};
        x_in = vertcat(u/100 + 0.5, theta/pi + 0.5);
        f_expl(3:4) = f_expl(3:4) + 0.01*bspline(x_in, C, knots, [3, 2], 2);
    end

    model = AcadosModel();
    model.f_impl_expr = xdot - f_expl;
    model.f_expl_expr = f_expl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.p = p;
    model.p_global = p_global;
    model.name = model_name;

end

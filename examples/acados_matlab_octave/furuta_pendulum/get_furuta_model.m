function model = get_furuta_model()
    import casadi.*

    % Distances
    L1 = 0.1035; % 103.5mm
    l2 = 0.0955; % 92.1mm

    % mass / inertia / gravity
    m2 = 0.192;        % kg
    J2 = 7.653e-04;    % kg*m^2
    g  = 9.81;         % m/s^2

    % inertia arm 1
    J1_ges = 5.3875e-04 + 0.75e-04; % J1 + m1*l1^2

    % inertia arm 2
    J2_ges = J2 + m2*l2^2;

    % total inertia at motor
    J0 = J1_ges + m2*L1^2;

    % damping hub motor
    b1 = 40*1e-4;

    % damping coupling between both arms
    k  = 0.098;
    b2 = 2*k*J2_ges;

    % applied torques
    tau2 = 0;

    % symbolic variables
    theta1  = SX.sym('theta1');
    theta2  = SX.sym('theta2');
    dtheta1 = SX.sym('dtheta1');
    dtheta2 = SX.sym('dtheta2');
    dtheta  = [dtheta1; dtheta2];
    tau1    = SX.sym('tau1');

    x    = [theta1; theta2; dtheta1; dtheta2];
    xdot = SX.sym('xdot', size(x));
    u    = tau1;

    % adjust coordinate (same as python: theta2 = theta2 - pi)
    theta2 = theta2 - pi;

    % trig helpers
    sin_theta_2    = sin(theta2);
    cos_theta_2    = cos(theta2);
    sin_2_theta_2  = sin(2*theta2);

    factor = m2*L1*l2;

    Matrix1 = [ J0 + J2_ges*sin_theta_2^2,    factor*cos_theta_2;
                factor*cos_theta_2,            J2_ges ];

    Matrix2 = [ b1 + 0.5*dtheta2*J2_ges*sin_2_theta_2,  0.5*dtheta1*J2_ges*sin_2_theta_2 - factor*sin_theta_2*dtheta2;
            -0.5*dtheta1*J2_ges*sin_2_theta_2,      b2 ];

    rhs = [tau1; tau2] - Matrix2 * [dtheta1; dtheta2] - [0; g*m2*l2*sin_theta_2];

    f_expl_expr = [dtheta; Matrix1 \ rhs];
    f_impl_expr = [ dtheta - xdot(1:2);
                    Matrix1*xdot(3:4) - rhs ];

    model = AcadosModel();
    model.name = 'furuta_model';
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.f_impl_expr = f_impl_expr;
    model.f_expl_expr = f_expl_expr;
end

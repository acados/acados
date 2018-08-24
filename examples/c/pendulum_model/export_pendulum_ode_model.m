function [ model ] = export_pendulum_ode_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr, f_expl_expr
    % that depends on the symbolic CasADi variables x, xdot, u, z,
    % and a model name, which will be used as a prefix for generated C
    % functions later on;
%% this model is based on the explicit pendulum model
% but formulated implicitly to test implicit integrators with it.
    %% CasADi
    import casadi.*
    if CasadiMeta.version()=='3.4.0'
        % casadi 3.4
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else
        % old casadi versions
        error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
    end
    model_name = 'pendulum_ode';
    

    %% Constants
    M = 1;
    m = 0.1;
    g = 9.81;
    l = 0.8;

    %% Set up States & Controls
    x1      = SX.sym('x1');
    theta   = SX.sym('theta');
    v1      = SX.sym('v1');
    dtheta  = SX.sym('dtheta');
    
    x = [x1; v1; theta; dtheta];

    % Controls
    F = SX.sym('F');
    u = F;
    
    %% xdot
    x1_dot      = SX.sym('x1_dot');
    theta_dot   = SX.sym('theta_dot');
    v1_dot      = SX.sym('v1_dot');
    dtheta_dot  = SX.sym('dtheta_dot');

    xdot = [ x1_dot; theta_dot; v1_dot; dtheta_dot ];
    
    %% algebraic variables
    z = [];
    
    %% Dynamics     
    denominator = M + m - m*cos(theta)*cos(theta);
    f_expl = [  v1; ...
                (-m*l*sin(theta)*dtheta*dtheta + m*g*cos(theta)*sin(theta)+F)/denominator; ...
                dtheta; ...
                (-m*l*cos(theta)*sin(theta)*dtheta*dtheta + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator)];
    
    f_impl = xdot - f_expl;
    
    model.f_impl_expr = f_impl;
    model.f_expl_expr = f_expl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.name = model_name;

end


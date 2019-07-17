function [ model ] = export_pure_linear_test_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr
    % that depends on the symbolic CasADi variables x, xdot, u, z,
    % and a model name, which will be used as a prefix for generated C
    % functions later on;
    
    %% CasADi
    import casadi.*
    casadi_version = CasadiMeta.version();
    if strcmp(casadi_version(1:3),'3.4') % require casadi 3.4.x
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else % old casadi versions
        error('Please download and install CasADi version 3.4.x to ensure compatibility with acados')
    end
    model_name_prefix = 'crane_nx9';
    

    %% Set up States & Controls
    xC = SX.sym('xC');     %States
    vC = SX.sym('vC');
    xL = SX.sym('xL');     
    vL = SX.sym('vL');
    uC = SX.sym('uC');
    uL = SX.sym('uL');
    theta = SX.sym('theta');
    omega = SX.sym('omega');
    q = SX.sym('q'); % a quadrature state
    x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega, q);
    xdot = SX.sym('xdot', size(x));

    uCR = SX.sym('uCR');  % Controls
    uLR = SX.sym('uLR');
    u = vertcat(uCR, uLR);

    z = SX.sym('z',0); % define an algebraic state;

    
    %% Dynamics: implicit DAE formulation (index-1)
    nx = length(x);
    nu = length(u);

    if 0
        A = rand(nx, nx);
        B = rand(nx, nu);
    else
        load('AB_test.mat')
    end
    f_expl = A*x + B*u;

    f_impl = (f_expl - xdot);

    model.f_impl_expr = f_impl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.name = model_name_prefix;
    
end
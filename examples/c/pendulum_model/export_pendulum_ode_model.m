function [ model ] = export_pendulum_ode_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr
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
    model_name_prefix = 'pendulum_ode_';
    

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
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.name = model_name_prefix;
end
% Sx = SX.sym('Sx',nx,nx);
% Sp = SX.sym('Sp',nx,nu);
% lambdaX = SX.sym('lambdaX',nx,1);
% 
% vdeX = SX.zeros(nx,nx);
% vdeX = vdeX + jtimes(f_expl,x,Sx);
% 
% vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u);
% vdeP = vdeP + jtimes(f_expl,x,Sp);
% 
% vdeFun = Function('vdeFun',{x,Sx,Sp,u},{f_expl,vdeX,vdeP});
% 
% jacX = SX.zeros(nx,nx) + jacobian(f_expl,x);
% jacFun = Function('jacFun',{x,u},{f_expl,jacX});
% 
% adj = jtimes(f_expl,[x;u],lambdaX,true);
% 
% adjFun = Function('adjFun',{x,lambdaX,u},{adj});
% 
% S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
% hess = S_forw.'*jtimes(adj,[x;u],S_forw);
% hess2 = [];
% for j = 1:nx+nu
%     for i = j:nx+nu
%         hess2 = [hess2; hess(i,j)];
%     end
% end
% 
% hessFun = Function('adjFun',{x,Sx,Sp,lambdaX,u},{adj,hess2});
% 
% opts = struct('mex', false, 'with_header', true, 'with_export', false, 'casadi_int', 'int');
% vdeFun.generate(['vde_forw_pendulum'], opts);
% jacFun.generate(['jac_pendulum'], opts);
% adjFun.generate(['vde_adj_pendulum'], opts);
% hessFun.generate(['vde_hess_pendulum'], opts);
% 
% p = vertcat(x1-l*sin(theta) - l, l*cos(theta) - l);
% quad_constraint = Function('position', {x}, {p, jacobian(p, x).T});
% quad_constraint.generate('position', opts);
% 
% h = mtimes(p.T, p);
% constraint = Function('constraint', {x}, {h, jacobian(h, x).T});
% constraint.generate('constraint', opts);

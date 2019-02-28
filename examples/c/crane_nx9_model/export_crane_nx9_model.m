function [ model ] = export_crane_nx9_model()
    %% this function generates an implicit ODE / index-1 DAE model,
    % which consists of a CasADi expression f_impl_expr
    % that depends on the symbolic CasADi variables x, xdot, u, z,
    % and a model name, which will be used as a prefix for generated C
    % functions later on;
    
    %% CasADi
    import casadi.*
    if CasadiMeta.version()=='3.4.0'
        % casadi 3.4
        casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
    else
        % old casadi versions
        error('Please download and install Casadi 3.4.0 to ensure compatibility with acados')
    end
    model_name_prefix = 'crane_nx9';

    %% Parameters (taken from Riens ACADO model)
    tau1 = 0.012790605943772;   a1   = 0.047418203070092;
    tau2 = 0.024695192379264;   a2   = 0.034087337273386;
    g = 9.81;

    %% Set up States & Controls
    xC = SX.sym('xC');     %States  
    vC = SX.sym('vC');
    xL = SX.sym('xL');     
    vL = SX.sym('vL');
    uC = SX.sym('uC');
    uL = SX.sym('uL');
    theta = SX.sym('theta');
    omega = SX.sym('omega'); 
    uCR = SX.sym('uCR');  % Controls
    uLR = SX.sym('uLR');
    q = SX.sym('q'); % a quadrature state

    %% explicit ODE formulation
    f_expl = vertcat(vC, ...
                      - 1/tau1 * (vC - a1 * uC), ...
                      vL,...
                      - 1/tau2 * (vL - a2 * uL), ...
                      uCR,...
                      uLR,...
                      omega, ...
                      - (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL, ...
                      uCR^2 + xL^2); % dynamics of quadrature state x2;

    %% Generalized nonlinear static feedback formulation (GNSF)
    x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega, q);
    u = vertcat(uCR, uLR);

    x1 = vertcat(xC, vC, xL, vL, uC, uL, theta, omega);
    x2 = q; % Linear output/quadrature state

    z = SX.sym('z',0,0); % define an algebraic state;
    xdot = SX.sym('xdot',length(x),1);
    % here z is just added to the ODE to have some algebraic states
    % z = (theta^2)/8
    f_impl = (f_expl - xdot);

    model.f_impl_expr = f_impl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.name = model_name_prefix;

end
% check_phi_jac_y_uhat = external('check_phi_jac_y_uhat', './check_phi_jac_y_uhat.so');

% 
% % get matrices - for use in final version
% model_matrices = SX.zeros(size([A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:)])) + ...
%     [A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:)];
% get_matrices_fun = Function('get_matrices_fun', {dummy}, {model_matrices(:)});
% get_matrices_fun.generate('get_matrices_fun', casadi_opts);
% 
% % generate Phi, f_LO
% f_lo_fun_jac_x1k1uz.generate(['f_lo_fun_jac_x1k1uz'], casadi_opts);
% Phi_inc_dy_fun.generate(['Phi_inc_dy_fun'], casadi_opts);

% Phi_fun.generate('Phi_fun',opts);
% f_fun.generate('f_fun', opts);
% 
% jac_Phi_y_fun.generate('jac_Phi_y_fun', opts);
% jac_Phi_u_fun.generate('jac_Phi_u_fun', opts);
% 
% jac_f_x1_fun.generate('jac_f_x1', opts);
% jac_f_u_fun.generate('jac_f_u',opts);
% 
% A_fun.generate('A_fun',opts);
% B_fun.generate('B_fun',opts);
% C_fun.generate('C_fun',opts);
% D_fun.generate('D_fun',opts);
% E_fun.generate('E_fun',opts);
% F_fun.generate('F_fun',opts);

% disp(['time to setup structure defining fcns/matrices and necessary derivatives = ', num2str(toc)]);

% save('structured_crane_model.mat','s')
% s.Phi_fun(ones(4,1), ones(2,1))
% % ode_fun = Function('ode_fun',{x,u},{f_expl});
% x_dot = SX.sym('x_dot',nx,1);
% f_impl = SX.zeros(nx,1)+(x_dot - f_expl);
% impl_odeFun = Function('impl_odeFun',{x_dot,x,u},{f_impl});
% 
% Sx = SX.sym('Sx',nx,nx);
% Sp = SX.sym('Sp',nx,nu);
% lambdaX = SX.sym('lambdaX',nx,1);
% 
% Derive Variational Differential Equations
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
% oj: The jtimes function optionally calculates the
% transposed-Jacobian-times-vector product, i.e.  reverse mode AD
% adj = jtimes(f_expl,[x;u],lambdaX,true);
% adj = jtimes(f_expl,[u;x],lambdaX,true);
% 
% adjFun = Function('adjFun',{x,lambdaX,u},{adj});
% 
% S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
% hess = S_forw.'*jtimes(adj,[x;u],S_forw);
% hess2 = [];
% for j = 1:nx+nu
%      for i = j:nx+nu
%          hess2 = [hess2; hess(i,j)];
%      end
% end
%  
% hessFun = Function('hessFun',{x,Sx,Sp,lambdaX,u},{adj,hess2});
% 
% opts = struct('mex', false);
% odeFun.generate(['ode_model'], opts);
% vdeFun.generate(['vde_forw_model'], opts);
% jacFun.generate(['jac_model'], opts);
% adjFun.generate(['vde_adj_model'], opts);
% hessFun.generate(['vde_hess_model'], opts);
% 
% implicit fcn generation for impl integrators
% x_dot = SX.sym('x_dot',nx,1);         
% f_impl = SX.zeros(nx,1)+(x_dot - f_expl);
% 
% impl_odeFun = Function('impl_odeFun',{x,x_dot,u},{f_impl});
% jac_x = SX.zeros(nx,nx) + jacobian(f_impl,x);
% % jac_xdot = SX.zeros(nx,nx) + jacobian(f_impl,x_dot);
% jac_u = SX.zeros(nx,nu) + jacobian(f_impl,u);
% % 
% impl_jacFun_x = Function('impl_jacFun_x',{x_dot,x,u},{jac_x});
% % impl_jacFun_xdot = Function('impl_jacFun_xdot',{x,x_dot,u},{jac_xdot});
% impl_jacFun_u = Function('impl_jacFun_u',{x_dot,x,u},{jac_u});

% 
% impl_odeFun.generate(['impl_ode'],opts);
% impl_jacFun_x.generate(['impl_jac_x'],opts);
% impl_jacFun_xdot.generate(['impl_jac_xdot'],opts);
% impl_jacFun_u.generate(['impl_jac_u'],opts);
% 
% x0 = zeros(nx,1); x0(3)=0.8;
% u0 = [40.108149413030752; -50.446662212534974];
% y0 = L_x *x0(1:nx1) + L_u * u0;
% k0 = 0*ones(nx,1);
% odeFun(x0,u0)

% full(impl_odeFun(x0,k0,u0))
% full(impl_jacFun_x(x0,k0,u0))
% full(impl_jacFun_xdot(x0,k0,u0))
% full(impl_jacFun_u(x0,k0,u0))
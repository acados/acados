%% export simple wind turbine model 
% model is taken from the paper:
% Real-time Economic Nonlinear Model Predictive Control for Wind Turbine
% Control - eq (1)

clc;
clear all;
close all;

% addpath('../../external/casadi-octave-v3.2.2')
import casadi.*


%% Parameters
Param = TurbineParameters;

rho = Param.General.rho;
A = pi*Param.Turbine.R^2;
r = Param.Turbine.i;
M = Param.Turbine.m_Te;
K = Param.Turbine.k_Te;
xi = Param.Turbine.c_Te;
J = Param.Turbine.J;
R = Param.Turbine.R;

w = 10; % OnlineData w; w0 = 10; sim_in.od = w0;

%% Set up States & Controls
Omega = SX.sym('Omega');     %States  
theta = SX.sym('theta');
dtheta= SX.sym('dtheta');     
T     = SX.sym('T');
x_w   = SX.sym('x_w');
x_wdot  = SX.sym('xdot');
uT    = SX.sym('uT');  %Controls
utheta = SX.sym('utheta');

% intermediate States
V = w - x_wdot;
lambda = R * Omega/V;
t2 = lambda.^2;
t3 = theta.^2;
Cp_poly = lambda.*5.138106966972632e-1 - t2.*5.364329820732732e-2 +...
    t3.*3.451233947321301e-3 + theta.*5.638964435900082e-2 + ...
    lambda.*t2.*1.768170297611305e-3 - lambda.*t3.*1.120255750841836e-3 -...
    lambda.*theta.*1.651646521371795e-2 + t2.*theta.*1.202000995895707e-3 +...
    t3.*theta.*6.300142515103073e-5 - 1.087868653202844;
Ct_poly = lambda.*4.730413934399427e-1 - t2.*4.333623415157917e-2 - t3.*2.240505670718569e-4 + ...
    theta.*7.820155434909468e-2 + lambda.*t2.*1.459194800744952e-3 - ...
    lambda.*t3.*2.559415177377994e-4 - lambda.*theta.*2.132624029767644e-2 + ...
    t2.*theta.*5.666961824200767e-4 + t3.*theta.*4.670564940602716e-5 - 9.343955026755452e-1;
TA = 0.5/J*rho*A*Cp_poly*V^3/Omega;
FA = 0.5/M*rho*A*Ct_poly*V^2;

%% Set up ODE
f_expl = vertcat( TA-T/(J*r), ...
             dtheta, ...
             utheta, ...
             uT, ...
             x_wdot, ...
             1/M*(-K * x_w -xi*x_wdot)+FA  );
         
x = vertcat( Omega, theta, dtheta, T, x_w, x_wdot);
u = vertcat( utheta, uT);

x1 = vertcat( Omega, theta, dtheta, T, x_w, x_wdot);
x2 = [];
z  = [];

nx1 = length(x1);
nu = length(u);
nx2 = length(x2);
nz = length(z);
nx = nx1 + nx2;

x1_dot = SX.sym('x1_dot',nx1,1);

opts = struct('mex', false);
     
%% structured
Phi = vertcat(TA, FA);
n_out = length(Phi);

A_mat = zeros(nx1+nz,nx1);
B = zeros(nx1+nz,nu );
C = zeros(nx1+nz, n_out);
E = eye(nx1+nz);

A_mat(1,4) = -1/(J*r); C(1,1) = 1;
A_mat(2,3) = 1;
B(3,1) = 1;
B(4,2) = 1;
A_mat(5,6) = 1;
A_mat(6,5) = -K/M; A_mat(6,6) = -xi/M; C(6,2)= 1;

f = SX.sym('f',0,0);

ALO = zeros(nx2);

% generate CasADi functions
jac_f_x1 = jacobian(f,x1);
jac_f_u = jacobian(f,u);
jac_f_z = jacobian(f,z);

Phi_fun = Function('Phi_fun', {x1_dot,x1,z,u}, {Phi});
f_fun = Function('f_los', {x1_dot,x1,z,u}, {f});

jac_f_x1_fun = Function('jac_f_x1_fun', {x1_dot,x1,z,u},{jac_f_x1});
jac_f_u_fun = Function('jac_f_u_fun', {x1_dot,x1,z,u},{jac_f_u});
jac_f_z_fun = Function('jac_f_z_fun', {x1_dot,x1,z,u},{jac_f_z});

s = struct('A', A_mat, 'B', B, 'C', C, 'E', E, 'ALO',ALO,...
    'Phi_fun', Phi_fun, 'f_fun', f_fun, ...
    'nx1', nx1, 'nx2', nx2, 'nu', nu, 'n_out', n_out, 'nx', nx, 'nz', nz,...
    'jac_f_x1_fun', jac_f_x1_fun, 'jac_f_u_fun', jac_f_u_fun, 'jac_f_z_fun', jac_f_z_fun);

EXPORT_ACADOS = 1;
if EXPORT_ACADOS
    odeFun = Function('odeFun',{x,u},{f_expl});
    x_dot = SX.sym('x_dot',nx,1);         
    f_impl = SX.zeros(nx,1)+(x_dot - f_expl); %% add SX.zeros to densify the output

    impl_odeFun = Function('impl_odeFun',{x,x_dot,u},{f_impl});
    jac_x = SX.zeros(nx,nx) + jacobian(f_impl,x);
    jac_xdot = SX.zeros(nx,nx) + jacobian(f_impl,x_dot);
    jac_u = SX.zeros(nx,nu) + jacobian(f_impl,u);
    
    Sx = SX.sym('Sx',nx,nx);
    Sp = SX.sym('Sp',nx,nu);
    lambdaX = SX.sym('lambdaX',nx,1);
    
    vdeX = SX.zeros(nx,nx);
    vdeX = vdeX + jtimes(f_expl,x,Sx);
    
    vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u);
    vdeP = vdeP + jtimes(f_expl,x,Sp);
    
    adj = jtimes(f_expl,[x;u],lambdaX,true);
    % adj = jtimes(f_expl,[u;x],lambdaX,true);

    adjFun = Function('adjFun',{x,lambdaX,u},{adj});

    S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
    hess = S_forw.'*jtimes(adj,[x;u],S_forw);
    hess2 = [];
    for j = 1:nx+nu
        for i = j:nx+nu
            hess2 = [hess2; hess(i,j)];
        end
    end

    hessFun = Function('hessFun',{x,Sx,Sp,lambdaX,u},{adj,hess2});
    
    vdeFun = Function('vdeFun',{x,Sx,Sp,u},{f_expl,vdeX,vdeP});
    
    jacX = SX.zeros(nx,nx) + jacobian(f_expl,x);
    jacFun = Function('jacFun',{x,u},{f_expl,jacX});
    
    odeFun.generate(['ode_model'], opts);
    
    vdeFun.generate(['vde_forw_model'], opts);
    jacFun.generate(['jac_model'], opts);
    adjFun.generate(['vde_adj_model'], opts);
    hessFun.generate(['vde_hess_model'], opts);

    impl_jacFun_x = Function('impl_jacFun_x',{x,x_dot,u},{jac_x});
    impl_jacFun_xdot = Function('impl_jacFun_xdot',{x,x_dot,u},{jac_xdot});
    impl_jacFun_u = Function('impl_jacFun_u',{x,x_dot,u},{jac_u});
    
    impl_odeFun.generate(['impl_ode'],opts);
    impl_jacFun_x.generate(['impl_jac_x'],opts);
    impl_jacFun_xdot.generate(['impl_jac_xdot'],opts);
    impl_jacFun_u.generate(['impl_jac_u'],opts);

    opts = struct('mex', false);
end
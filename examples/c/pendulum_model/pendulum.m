clc;
clear all;
close all;

GENERATE_LQR_GAIN = 1;

addpath('../../external/casadi-octave-v3.2.2')
import casadi.*

% variables
x1 = SX.sym('x1');
theta = SX.sym('theta');
v1 = SX.sym('v1');
dtheta = SX.sym('dtheta');

F = SX.sym('F');

x = [x1; theta; v1; dtheta];
u = F;

nx = length(x);
nu = length(u);

% ODE system
f_expl = dyn(x, u);


         
Sx = SX.sym('Sx',nx,nx);
Sp = SX.sym('Sp',nx,nu);
lambdaX = SX.sym('lambdaX',nx,1);

vdeX = SX.zeros(nx,nx);
vdeX = vdeX + jtimes(f_expl,x,Sx);

vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u);
vdeP = vdeP + jtimes(f_expl,x,Sp);

vdeFun = Function('vdeFun',{x,Sx,Sp,u},{f_expl,vdeX,vdeP});

jacX = SX.zeros(nx,nx) + jacobian(f_expl,x);
jacFun = Function('jacFun',{x,u},{f_expl,jacX});

adj = jtimes(f_expl,[x;u],lambdaX,true);

adjFun = Function('adjFun',{x,lambdaX,u},{adj});

S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
hess = S_forw.'*jtimes(adj,[x;u],S_forw);
hess2 = [];
for j = 1:nx+nu
    for i = j:nx+nu
        hess2 = [hess2; hess(i,j)];
    end
end

hessFun = Function('adjFun',{x,Sx,Sp,lambdaX,u},{adj,hess2});

opts = struct('mex', false);
vdeFun.generate(['vde_forw_pendulum'], opts);
jacFun.generate(['jac_pendulum'], opts);
adjFun.generate(['vde_adj_pendulum'], opts);
hessFun.generate(['vde_hess_pendulum'], opts);

if GENERATE_LQR_GAIN % generate LQR gain
   
    Ts = 0.01;
    
    % fixed step Runge-Kutta 4 integrator
    M = 1; % RK4 steps per interval
    DT = Ts;
    X0 = SX.sym('X0', 4,1);
    U = SX.sym('U',1,1);
    X = X0;
    for j=1:M
        [k1] = dyn(X, U);
        [k2] = dyn(X + DT/2 * k1, U);
        [k3] = dyn(X + DT/2 * k2, U);
        [k4] = dyn(X + DT * k3, U);
        X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    F = Function('F', {X0, U}, {X});
        
    J_x = Function('J_x', {X0, U}, {jacobian(X, X0)});
    J_u = Function('J_u', {X0, U}, {jacobian(X, U)});
    
    A = full(J_x(zeros(4,1), zeros(1,1)));
    B = full(J_u(zeros(4,1), zeros(1,1)));
    
    Q = zeros(4,4);
    R = zeros(1,1);
    
    Q(1,1) = 1e-1;
    Q(2,2) = 1.0;
    Q(3,3) = 0.1;
    Q(4,4) = 2e-3;

    R(1,1) = 5e-4;
    
    [K, P] = dlqr(A, B, Q, R)
    
end

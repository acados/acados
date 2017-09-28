clc;
clear all;
close all;

GENERATE_LQR_GAIN = 0;

addpath('../../external/casadi-octave-v3.2.2')
import casadi.*

% variables
q = SX.sym('q', 4, 1);
omega = SX.sym('omega', 3, 1);
W = SX.sym('W', 4, 1);
rW = SX.sym('rW', 4, 1);

x = [q; omega; W];
u = rW;

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
vdeFun.generate(['vde_forw_quadcopter'], opts);
jacFun.generate(['jac_quadcopter'], opts);
adjFun.generate(['vde_adj_quadcopter'], opts);
hessFun.generate(['vde_hess_quadcopter'], opts);

if GENERATE_LQR_GAIN % generate LQR gain

    Ts = 0.01;

    % fixed step Runge-Kutta 4 integrator
    M = 1; % RK4 steps per interval
    DT = Ts;
    X0 = SX.sym('X0', 11,1);
    U = SX.sym('U',4,1);
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

    A = full(J_x(zeros(11,1), zeros(4,1)));
    B = full(J_u(zeros(11,1), zeros(4,1)));

    Q = eye(11);

    R = eye(4);

    [K, P] = dlqr(A, B, Q, R)

end

% generate code for residuals
nr = 15;
nr_end = 11;
res_exp = [x;u];
res_end_exp = [x];

jac_res_exp = SX.zeros(nr,nx+nu) + jacobian(res_exp,[x;u]);
ls_res_Fun = Function('ls_res_Fun', {x,u}, {res_exp,jac_res_exp});
jac_res_end_exp = SX.zeros(nr_end,nx) + jacobian(res_end_exp,[x]);
ls_res_end_Fun = Function('ls_res_end_Fun', {x,u}, {res_end_exp,jac_res_end_exp});

ls_res_Fun.generate(['ls_res_quadcopter'], opts);
ls_res_end_Fun.generate(['ls_res_end_quadcopter'], opts);

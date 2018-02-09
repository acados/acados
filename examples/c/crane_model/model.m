clc;
clear all;
close all;

import casadi.*

% constants
m = 0.1;
g = 9.81;
l = 0.8;

% variables
x1 = SX.sym('x1');
theta = SX.sym('theta');
v1 = SX.sym('v1');
dtheta = SX.sym('dtheta');

F = SX.sym('F');

M = SX.sym('M');

x = [x1; theta; v1; dtheta];
u = F;
p = M;

nx = length(x);
nu = length(u);
np = length(p);

% ODE system

f_expl = [   v1; ...
             dtheta; ...
             (- l*m*sin(theta)*dtheta^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta)^2); ...
             (- l*m*cos(theta)*sin(theta)*dtheta^2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta)^2)) ];
        
odeFun = Function('odeFun',{x,u,p},{f_expl});

Sx = SX.sym('Sx',nx,nx);
Sp = SX.sym('Sp',nx,nu);
lambdaX = SX.sym('lambdaX',nx,1);

vdeX = jtimes(f_expl,x,Sx);

vdeP = jacobian(f_expl,u);
vdeP = vdeP + jtimes(f_expl,x,Sp);

vdeFun = Function('vdeFun',{x,Sx,Sp,u,p},{f_expl,vdeX,vdeP});

jacX = jacobian(f_expl,x);
jacFun = Function('jacFun',{x,u,p},{f_expl,jacX});

adj = jtimes(f_expl,[x;u],lambdaX,true);
% adj = jtimes(f_expl,[u;x],lambdaX,true);

adjFun = Function('adjFun',{x,lambdaX,u,p},{adj});

S_forw = vertcat(horzcat(Sx, Sp), horzcat(zeros(nu,nx), eye(nu)));
hess = S_forw.'*jtimes(adj,[x;u],S_forw);
hess2 = [];
for j = 1:nx+nu
    for i = j:nx+nu
        hess2 = [hess2; hess(i,j)];
    end
end

hessFun = Function('hessFun',{x,Sx,Sp,lambdaX,u,p},{adj,hess2});

opts = struct('mex', false);
odeFun.generate(['ode_model'], opts);
vdeFun.generate(['vde_forw_model'], opts);
jacFun.generate(['jac_model'], opts);
adjFun.generate(['vde_adj_model'], opts);
hessFun.generate(['vde_hess_model'], opts);

x_dot = SX.sym('x_dot',nx,1);         
f_impl = (x_dot - f_expl); %% add SX.zeros to densify the output

impl_odeFun = Function('impl_odeFun',{x,x_dot,u,p},{f_impl});
jac_x = jacobian(f_impl,x);
jac_xdot = jacobian(f_impl,x_dot);
jac_u = jacobian(f_impl,u);

impl_jacFun = Function('impl_jacFun',{x,x_dot,u,p},{jac_x,jac_xdot,jac_u});

opts = struct('mex', false);
% 
impl_odeFun.generate(['impl_ode'],opts);
impl_jacFun.generate(['impl_jac'],opts);

x0 = [0;3.14;0;0];
u0 = 0;
k0 = 0*ones(4,1);

% odeFun(x0,u0)

impl_odeFun(x0,k0,u0,1)

impl_jacFun(x0,k0,u0,1)
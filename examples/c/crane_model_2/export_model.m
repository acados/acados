clc;
clear all;
close all;

import casadi.*

%% Parameters (taken from Riens ACADO model)
tau1 = 0.012790605943772;   a1   = 0.047418203070092;
tau2 = 0.024695192379264;   a2   = 0.034087337273386;
g = 9.81;

%% Set up States & Controls
xC = SX.sym('xC');  %States  
vC = SX.sym('vC');
xL = SX.sym('xL');
vL = SX.sym('vL');
uC = SX.sym('uC');
uL = SX.sym('uL');
theta = SX.sym('theta');
omega = SX.sym('omega'); 
uCR = SX.sym('uCR');  %Controls
uLR = SX.sym('uLR');

q = SX.sym('q');
x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega, q);
u = vertcat(uCR, uLR);

x1 = vertcat(xC, vC, xL, vL, uC, uL, theta, omega);
x2 = q; % some quadrature state

z = SX.sym('z'); % define an algebraic state
nz = length(z);

f_expl = vertcat(vC, ...
                  - 1/tau1 * (vC - a1 * uC), ...
                  vL,...
                  - 1/tau2 * (vL - a2 * uL), ...
                  uCR,...
                  uLR,...
                  omega, ...
                  - (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL, ...
                  uCR^2 + xL^2); % dynamics of quadrature state x2;

f = uCR^2 + xL^2;

nx1 = length(x1);
nu = length(u);
nx2 = length(x2);
nx = nx1 + nx2;

x1_dot = SX.sym('x1_dot',nx1,1);

% % constants
% M = 1;
% m = 0.1;
% g = 9.81;
% l = 0.8;
% 
% % variables
% x1 = SX.sym('x1');
% theta = SX.sym('theta');
% v1 = SX.sym('v1');
% dtheta = SX.sym('dtheta');
% 
% F = SX.sym('F');
% 
% x = [x1; theta; v1; dtheta];
% u = F;
% 
% nx = length(x);
% nu = length(u);
% 
% % ODE system
% 
% f_expl = [   v1; ...
%              dtheta; ...
%              (- l*m*sin(theta)*dtheta^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta)^2); ...
%              (- l*m*cos(theta)*sin(theta)*dtheta^2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta)^2)) ];


        
% odeFun = Function('odeFun',{x,u},{f_expl});

x_dot = SX.sym('x_dot',nx,1);         
f_impl = (x_dot - f_expl); %% add SX.zeros to densify the output

jac_x =  jacobian(f_impl,x);
jac_xdot =  jacobian(f_impl,x_dot);
jac_u = jacobian(f_impl,u);

impl_odeFun = Function('impl_odeFun',{x,x_dot,u},{f_impl});
impl_jacFun_x = Function('impl_jacFun_x',{x,x_dot,u},{jac_x});
impl_jacFun_xdot = Function('impl_jacFun_xdot',{x,x_dot,u},{jac_xdot});
impl_jacFun_u = Function('impl_jacFun_u',{x,x_dot,u},{jac_u});

opts = struct('mex', false);
% 
impl_odeFun.generate(['impl_ode'],opts);
impl_jacFun_x.generate(['impl_jac_x'],opts);
impl_jacFun_xdot.generate(['impl_jac_xdot'],opts);
impl_jacFun_u.generate(['impl_jac_u'],opts);
function [ model ] = TEST_export_inverted_pendulum_dae_model_MX()
%% this function generates an implicit ODE / index-1 DAE model,
% which consists of a CasADi expression f_impl_expr
% that depends on the symbolic CasADi variables x, xdot, u, z,
% and a model name, which will be used as a prefix for generated C
% functions later on;

model_name_prefix = 'inverted_pendulum_MX';

%% CasADi
import casadi.*
if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0')
end

%% Parameters (taken from Rien Quirynens Master Thesis)
m = 2;
g = 9.81;
M = 3.5;
I = 0.1;

%% Set up States & Controls
xpos    = MX.sym('xpos');     % Differential States
ypos    = MX.sym('ypos');
alpha   = MX.sym('alpha');     
vx      = MX.sym('vx');
vy      = MX.sym('vy');
valpha  = MX.sym('valpha');
x = vertcat(xpos, ypos, alpha, vx, vy, valpha);

ax      = MX.sym('ax');     % Algebraic states
ay      = MX.sym('ay');
aalpha  = MX.sym('aalpha');
Fx      = MX.sym('Fx');
Fy      = MX.sym('Fy');
z = vertcat(ax, ay, aalpha, Fx, Fy);

%% controls
u       = MX.sym('u');  % Controls

%% xdot
xpos_dot    = MX.sym('xpos_dot');     % Differential States
ypos_dot    = MX.sym('ypos_dot');
alpha_dot   = MX.sym('alpha_dot');     
vx_dot      = MX.sym('vx_dot');
vy_dot      = MX.sym('vy_dot');
valpha_dot  = MX.sym('valpha_dot');

xdot = [xpos_dot; ypos_dot; alpha_dot; vx_dot; vy_dot; valpha_dot];

%% parameters
p = MX.sym('parameters',0,0);

%% dimensions
nx = length(x);
nu = length(u);
nz = length(z);
np = length(p);

%% Dynamics: implicit DAE formulation (index-1)
% x = vertcat(xpos, ypos, alpha, vx, vy, valpha);
% z = vertcat(ax, ay, aalpha, Fx, Fy);
f_impl = vertcat(xpos_dot - vx, ...
                 ypos_dot - vy, ...
                 alpha_dot - valpha, ...
                 vx_dot - ax, ...
                 vy_dot - ay, ...
                 valpha_dot - aalpha, ...
                 m * ax - (Fx + u), ...
                 m * ay + m * g - Fy, ...
                 I * aalpha - M - (Fx + u) * ypos + Fy * xpos, ...
                 ax + vy * valpha + ypos * aalpha, ...
                 ay - vx * valpha - xpos * aalpha);

%                  ay - vx * vx - xpos * aalpha);

%% initial value
x0 = [1; -5; 1; 0.1; -0.5; 0.1];
z0 = [-1.5; -0.3; -0.3; -3; 19];
u0 = 1;


model.f_impl_expr = f_impl;
model.x = x;
model.xdot = xdot;
model.u = u;
model.z = z;
model.p = p;
model.name = model_name_prefix;


%% test model function -- check equivalence
CHECK = 0;
if CHECK
x0dot = zeros(nx,1);
x0dot = rand(nx,1);
% x0    = rand(nx,1);
z0    = rand(nz,1);

x1_0 = x0(1:nx1);
x1dot_0 = x0dot(1:nx1);
x2_0 = x0(nx1+1:nx);

y0 = L_x * x1_0 + L_xdot * x1dot_0 + L_z * z0;

uhat0 = L_u * u0;
phi_val0 = phi_fun(y0, uhat0);
gnsf_val0 = [E * [x1dot_0; z0] - full(A*x1_0 + B * u0 + C * phi_val0 + c); 
            full(f_lo_fun(x1_0, x1dot_0, z0, u0))]

f_imp_0 = full(impl_ode_fun(x0, x0dot, u0, z0))
gnsf_val0 - f_imp_0

assert(norm(gnsf_val0 - f_imp_0) < 1e-14)
end






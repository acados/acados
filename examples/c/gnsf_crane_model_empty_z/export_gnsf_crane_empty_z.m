clc;
clear;
close all;

% addpath('../../external/casadi-octave-v3.2.2')
import casadi.*

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
uCR = SX.sym('uCR');  %Controls
uLR = SX.sym('uLR');

q = SX.sym('q');
x = vertcat(xC, vC, xL, vL, uC, uL, theta, omega, q);
u = vertcat(uCR, uLR);

x1 = vertcat(xC, vC, xL, vL, uC, uL, theta, omega);
x2 = q; % some quadrature state

z = SX.sym('z',0); % define an algebraic state
size_z = size(z);
if min(size_z) == 0
    nz = 0;
else nz = max(size_z);
end

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

casadi_opts = struct('mex', false);

%% structured
tic
A = zeros(nx1+nz,nx1);
A(1,2) = 1;
A(2,2) = -1/tau1;
A(2,5) = a1/tau1;
A(3,4) = 1;
A(4,4) = -1/tau2;
A(4,6) = a2/tau2;
A(7,8) = 1;

B = zeros( nx1+nz, nu);
B(5,1) = 1;
B(6,2) = 1;

E = eye(nx1+nz);

Phi = [- (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL];%- (a1 * uCR * cos(y(3)) + g* sin(y(3)) + y(2)*y(4) ) / y(1);

n_out = length(Phi);

C = zeros( nx1+nz, n_out); C(8,1) = 1; %C(9,2) = 1;

ALO = zeros(nx2);

% f = uCR^2 + xL^2;
jac_f_x1 = jacobian(f,x1);
jac_f_k1 = jacobian(f,x1_dot);
jac_f_u = jacobian(f,u);
jac_f_z = jacobian(f,z);

Phi_fun = Function('Phi_fun', {x1_dot,x1,z,u}, {Phi});
f_fun = Function('f_los', {x1_dot,x1,z,u}, {f});

s = struct('A', A, 'B', B, 'C', C, 'E', E, 'ALO',ALO,...
    'Phi_fun', Phi_fun, 'f_fun', f_fun, ...
    'nx1', nx1, 'nx2', nx2, 'nu', nu, 'n_out', n_out, 'nx', nx, 'nz', nz);
% save('structured_crane_model.mat','s')
% s.Phi_fun(ones(4,1), ones(2,1))
% ode_fun = Function('ode_fun',{x,u},{f_expl});
x_dot = SX.sym('x_dot',nx,1);

EXPORT_ACADOS = 1;

%% choose options
q = 4; % number of stages;
stepsize = .1;
n_steps = 2;
disp(['Number of stages of RK method = ' num2str(q)]);
disp(['Number of steps per simulation =   ' num2str(n_steps)]);
adj_sens_mode = 1;
forw_sens_mode = 1;
opts = set_integrator_opts(q, stepsize, n_steps, forw_sens_mode, adj_sens_mode);
opts.sensi_mode = 'direct';
opts.max_newton = 10;


%% Using the GNLF_integrator
I_stages = eye( opts.n_stages );
%% generate submatrices
E11 = s.E(1:nx1, 1:nx1);
E12 = s.E(1:nx1, 1+nx1:nx1+nz);
E21 = s.E(1+nx1:nx1+nz, 1:nx1);
E22 = s.E(1+nx1:nx1+nz, 1+nx1:nx1+nz);

A1 = s.A(1:nx1, :); A2 = s.A(nx1+1:nx1+nz, :);
B1 = s.B(1:nx1, :); B2 = s.B(nx1+1:nx1+nz, :);
C1 = s.C(1:nx1, :); C2 = s.C(nx1+1:nx1+nz, :);

% generate fat matrices
AA1 = repmat(A1,opts.n_stages,1);
AA2 = repmat(A2,opts.n_stages,1);
BB1 = repmat(B1,opts.n_stages,1);
BB2 = repmat(B2,opts.n_stages,1);
CC1 = kron(I_stages,C1);
CC2 = kron(I_stages,C2);
DD1 = -kron(I_stages, E12);
DD2 = opts.dt * kron(opts.A_butcher, A2) - kron(I_stages, E21);
EE1 = kron(I_stages, E11) - opts.dt * kron(opts.A_butcher, A1);
EE2 = kron(I_stages, E22);

% needed to get K-values from f,u,x
EE1inv = inv(EE1);
EE2inv = inv(EE2);

PP1 = inv(eye(nx1*opts.n_stages) - EE1inv * DD1 * EE2inv  * DD2); 
PP2 = EE1inv * DD1 * EE2inv;
KKf = PP1 * ( PP2 * CC2 + EE1inv * CC1);
KKu = PP1 * ( PP2 * BB2 + EE1inv * BB1);
KKx = PP1 * ( PP2 * AA2 + EE1inv * AA1);

% to get Z-values
ZZf = EE2inv * (DD2 * KKf + CC2);
ZZu = EE2inv * (DD2 * KKu + BB2);
ZZx = EE2inv * (DD2 * KKx + AA2);

% needed for linear output system:
M2 = eye( s.nx2 * opts.n_stages);
M2 =  M2 - kron( opts.dt * opts.A_butcher, s.ALO);
M2inv = M2^-1;

% put precomputed matrices in struct
s.KKf = KKf;
s.KKu = KKu;
s.KKx = KKx;

s.ZZf = ZZf;
s.ZZu = ZZu;
s.ZZx = ZZx;
%% generate residual function
ff = SX.sym('ff', n_out * opts.n_stages,1);
x1_init = SX.sym('x1_init', nx1, 1);
u_in = SX.sym('u_in', nu, 1);
%     keyboard
K1 = KKf * ff + KKu * u_in + KKx * x1_init;
Z  = ZZf * ff + ZZu * u_in + ZZx * x1_init;
stageval = kron(opts.A_dt, eye(nx1)) * K1 + repmat( x1_init, opts.n_stages, 1);
res = SX.sym('res', opts.n_stages * n_out ,1);
for ii = 1 : opts.n_stages
    ind_out = index( n_out, ii);
    ind_x1  = index( nx1,   ii);
    ind_z   = index( nz ,   ii);
    res( ind_out,1 ) = ff(ind_out) - s.Phi_fun( K1( ind_x1 ), stageval(ind_x1), Z(ind_z), u_in);
end
jac_res_ff = jacobian(res, ff);


res_inc_Jff = [res, jac_res_ff]; % SX.zeros( q*n_out, q*n_out +1) +
% stageval_fun = Function('stageval_fun', {ff, x1_init, u_in}, {stageval}); %% TODO: maybe use

res_inc_Jff_fun = Function('res_inc_Jff_fun', {ff, x1_init, u_in}, {res_inc_Jff});
% jac_res_ffx1u = jacobian(res, [ff; x1_init; u_in]);
% this is more efficient, because of backward AD
jac_res_ffx1u = transpose(jtimes(res, [ff; x1_init; u_in], eye(opts.n_stages * n_out), true)); %SX.zeros(size(jac_res_ffx1u)) +
jac_res_ffx1u_fun = Function('jac_res_ffx1u_fun',{ff, x1_init, u_in}, {jac_res_ffx1u});

res_inc_Jff_fun.generate(['res_inc_Jff_fun'], casadi_opts);
jac_res_ffx1u_fun.generate(['jac_res_ffx1u_fun'],casadi_opts);



%% generate f_LO function
f_LO_inc_J_x1k1uz =  [f, jac_f_x1, jac_f_k1, jac_f_z, jac_f_u]; % + SX.zeros(nz, 1 + 2*nx1 + nz + nu)
f_LO_inc_J_x1k1uz_fun = Function('f_LO_inc_J_x1k1uz_fun', {x1, x1_dot, z, u}, ...
    {f_LO_inc_J_x1k1uz});
f_LO_inc_J_x1k1uz_fun.generate(['f_LO_inc_J_x1k1uz_fun'], casadi_opts);

%% generate matrices as functions:
dummy = SX.sym('dummy');
KKmatrices = SX.zeros( nx1 * q, n_out * q + nu + nx1 );
KKmatrices = KKmatrices + [KKf, KKx, KKu];

ZZmatrices = SX.zeros( nz * q, n_out * q + nu + nx1 );
ZZmatrices = ZZmatrices + [ZZf, ZZx, ZZu];

Butcher = SX.zeros( q, q+2) + [ opts.A_dt, opts.b_dt, opts.c_butcher];

dK2_dx2 = M2inv * repmat(s.ALO,opts.n_stages,1);
ALO_M2_dK2dx2 = SX.zeros(size([ ALO(:); M2inv(:); dK2_dx2(:)])) + [ALO(:); M2inv(:); dK2_dx2(:) ];

But_KK_ZZ_LO = [ Butcher(:); KKmatrices(:); ZZmatrices(:); ALO_M2_dK2dx2];
But_KK_ZZ_LO_fun = Function('But_KK_ZZ_LO_fun',{dummy}, {But_KK_ZZ_LO});
But_KK_ZZ_LO_fun.generate('But_KK_ZZ_LO_fun', casadi_opts);

%% generate functions
ints = SX.zeros(8,1) + [s.nx, s.nu, s.nz, s.nx1, s.nx2, q, n_steps, s.n_out]';
get_ints_fun = Function('get_ints_fun',{x},{[s.nx, s.nu, s.nz, s.nx1, s.nx2, q, n_steps, s.n_out]});
    get_ints_fun.generate('get_ints_fun', casadi_opts);

ff0 = zeros( s.n_out * q,1);
x1_0 = [0 0 0.8 0 0 0 0 0]';
u0 = zeros( s.nu, 1);
res_inc_Jff_val = full(res_inc_Jff_fun( ff0, x1_0, u0));

disp(['time in export_gnsf_crane = ',num2str(toc)]); 

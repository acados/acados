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

opts = struct('mex', false);

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

%% choose options
q = 4; % number of stages;
stepsize = .1;
n_steps = 2;
disp(['Number of stages of RK method = ' num2str(q)]);
disp(['Number of steps per simulation =   ' num2str(n_steps)]);
adj_sens_mode = 1;
forw_sens_mode = 1;
opts = set_integrator_opts(q, stepsize, n_steps, forw_sens_mode, adj_sens_mode);


%% nonlinearity function
Phi_old = [- (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL;
      (theta^2)/8];%- (a1 * uCR * cos(y(3)) + g* sin(y(3)) + y(2)*y(4) ) / y(1);
y = vertcat(xL, vL, theta, omega, uCR);

n_in = length(y);

% L_xdot = zeros(n_in,nx1); % do jacobian(y,x_dot);
% L_x    = zeros(n_in,nx1); L_x(1,3) = 1; L_x(2,4) = 1; L_x(3,7) = 1; L_x(4,8) = 1;
% L_z    = zeros(n_in, nz);
% L_u    = zeros(n_in, nu); L_u(5,1) = 1;

L_x_fun = Function('L_x_fun',{x1},{jacobian(y,x1)});
L_xdot_fun = Function('L_x_fun',{x1},{jacobian(y,x1_dot)});
L_z_fun = Function('L_z_fun',{x1},{jacobian(y,z)});
L_u_fun = Function('L_u_fun',{x1},{jacobian(y,u)});

L_x = full(L_x_fun(0));
L_xdot = full(L_xdot_fun(0));
L_u = full(L_u_fun(0));
L_z = full(L_z_fun(0));

y_check = L_xdot * x1_dot +L_x * x1 + L_z * z + L_u * u; %% THis should be the same as y

Phi = Phi_old;
n_out = length(Phi);

C = zeros( nx1+nz, n_out); C(8,1) = 1; C(9,2) = 1;

Phi_fun = Function('Phi_fun',{y}, {Phi});

jac_Phi_y = jacobian(Phi,y);
Phi_inc_dy_fun = Function('Phi_inc_dy_fun', {y}, {[Phi, jac_Phi_y]});
jac_Phi_y_fun = Function('jac_Phi_y_fun', {y}, {[jac_Phi_y]});

%% Linear output
ALO = zeros(nx2);
% A2(1,1) = 1;

% f = uCR^2 + xL^2;
jac_f_x1 = jacobian(f,x1);
jac_f_u = jacobian(f,u);
jac_f_z = jacobian(f,z);
jac_f_k1 = jacobian(f,x1_dot);

f_fun = Function('f_los', {x1_dot,x1,z,u}, {f});

% jac_Phi_u_fun = Function('jac_Phi_u_fun', {y,u},{jac_Phi_u});
% 
jac_f_x1_fun = Function('jac_f_x1_fun', {x1_dot,x1,z,u},{jac_f_x1});
jac_f_u_fun = Function('jac_f_u_fun', {x1_dot,x1,z,u},{jac_f_u});
jac_f_z_fun = Function('jac_f_z_fun', {x1_dot,x1,z,u},{jac_f_z});

f_LO_inc_J_x1k1uz =  [f, jac_f_x1, jac_f_k1, jac_f_z, jac_f_u]; % + SX.zeros(nz, 1 + 2*nx1 + nz + nu)
f_LO_inc_J_x1k1uz_fun = Function('f_LO_inc_J_x1k1uz_fun', {x1, x1_dot, z, u}, ...
    {f_LO_inc_J_x1k1uz});

% A_fun = Function('A_fun', {x1}, {A});
s = struct('A', A, 'B', B, 'C', C, 'E', E, 'ALO',ALO, 'L_x', L_x, 'L_xdot', L_xdot, 'L_z', L_z, 'L_u', L_u, ...
    'Phi_inc_dy_fun', Phi_inc_dy_fun, 'jac_Phi_y_fun', jac_Phi_y_fun, 'f_fun', f_fun, 'n_in', n_in,...
    'nx1', nx1, 'nx2', nx2, 'nu', nu, 'n_out', n_out, 'nx', nx, 'nz', nz,...
    'jac_f_x1_fun', jac_f_x1_fun, 'jac_f_u_fun', jac_f_u_fun, 'jac_f_z_fun', jac_f_z_fun, ...
    'f_LO_inc_J_x1k1uz_fun', f_LO_inc_J_x1k1uz_fun);

%% generate functions
casadi_opts = struct('mex', false);
dummy = SX.sym('dummy');

% get ints
get_ints_fun = Function('get_ints_fun',{x},{[s.nx, s.nu, s.nz, s.nx1, s.nx2, q, n_steps, s.n_out, s.n_in]});
get_ints_fun.generate('get_ints_fun', casadi_opts);

% get matrices - for use in final version
model_matrices = [A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:)];
get_matrices_fun = Function('get_matrices_fun', {x}, {model_matrices(:)});
get_matrices_fun.generate('get_matrices_fun', casadi_opts);

% generate Phi, f_LO
f_LO_inc_J_x1k1uz_fun.generate(['f_LO_inc_J_x1k1uz_fun'], casadi_opts);
Phi_inc_dy_fun.generate(['Phi_inc_dy_fun'], casadi_opts);

%% generate fat matrices
%generate submatrices
I_stages = eye( q );
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

PP1 = inv(eye(nx1*opts.n_stages) - (EE1 \ DD1) * EE2inv  * DD2); 
PP2 = (EE1 \ DD1) * EE2inv;
KKf = PP1 * ( PP2 * CC2 + EE1 \ CC1);
KKu = PP1 * ( PP2 * BB2 + EE1 \ BB1);
KKx = PP1 * ( PP2 * AA2 + EE1 \ AA1);

% to get Z-values
ZZf = EE2 \ (DD2 * KKf + CC2);
ZZu = EE2 \ (DD2 * KKu + BB2);
ZZx = EE2 \ (DD2 * KKx + AA2);

% to get yy-values:
LLZ = kron(I_stages, s.L_z);
LLx = repmat(s.L_x, opts.n_stages, 1);
LLK = kron(opts.dt * opts.A_butcher, s.L_x) + kron(I_stages, s.L_xdot);
LLu = repmat(s.L_u, opts.n_stages, 1);

YYx = LLx + LLK * KKx + LLZ * ZZx;
YYu = LLu + LLK * KKu + LLZ * ZZu;
YYf = LLK * KKf + LLZ * ZZf;

s.YYx = YYx;
s.YYu = YYu;
s.YYf = YYf;

% needed for linear output system:
M2 = eye( s.nx2 * opts.n_stages);
M2 =  M2 - kron( opts.dt * opts.A_butcher, s.ALO);
M2inv = M2^-1;

KKmatrices = SX.zeros( nx1 * q, n_out * q + nu + nx1 );
KKmatrices = KKmatrices + [KKf, KKx, KKu];

YYmatrices = SX.zeros( n_in*q, n_out*q + nx1 +nu) + [YYf, YYx, YYu];

ZZmatrices = SX.zeros( nz * q, n_out * q + nu + nx1 );
ZZmatrices = ZZmatrices + [ZZf, ZZx, ZZu];

Butcher = SX.zeros( q, q+2) + [ opts.A_dt, opts.b_dt, opts.c_butcher];

dK2_dx2 = M2inv * repmat(s.ALO,opts.n_stages,1);
ALO_M2_dK2dx2 = SX.zeros(size([ ALO(:); M2inv(:); dK2_dx2(:)])) + [ALO(:); M2inv(:); dK2_dx2(:) ];

But_KK_YY_ZZ_LO = [ Butcher(:); KKmatrices(:); YYmatrices(:); ZZmatrices(:); ALO_M2_dK2dx2];
But_KK_YY_ZZ_LO_fun = Function('But_KK_YY_ZZ_LO_fun',{dummy}, {But_KK_YY_ZZ_LO});
But_KK_YY_ZZ_LO_fun.generate('But_KK_YY_ZZ_LO_fun', casadi_opts);



disp(['time in export_gnsf_crane = ',num2str(toc)]); 


% Phi_fun.generate('Phi_fun',opts);
% f_fun.generate('f_fun', opts);

% save('structured_crane_model.mat','s')
% s.Phi_fun(ones(4,1), ones(2,1))
% ode_fun = Function('ode_fun',{x,u},{f_expl})
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

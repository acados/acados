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

z = SX.sym('z',2); % define an algebraic state;
% here z is just added to the ODE to have some algebraic states
% z = (theta^2)/8 + xL
nz = length(z);

f = uCR^2 + xL^2;

nx1 = length(x1);
nu = length(u);
nx2 = length(x2);
nx = nx1 + nx2;

x1_dot = SX.sym('x1_dot',nx1,1);

if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0')
end
model_name_prefix = 'crane_dae_';

%% Model defining matrices
A = zeros(nx1+nz,nx1);
A(1,2) = 1;
A(2,2) = -1/tau1;
A(2,5) = a1/tau1;
A(3,4) = 1;
A(4,4) = -1/tau2;
A(4,6) = a2/tau2;
A(7,8) = 1;

B = zeros( nx1 + nz, nu);
B(5,1) = 1;
B(6,2) = 1;

E = eye(nx1+nz);

c = zeros(nx1 + nz,1);
%% nonlinearity function
% gather all nonlinear parts
phi = [- (a1 * uCR * cos(theta) + g* sin(theta) + 2*vL*omega) / xL;
      (theta^2)/8 + xL;
      cos(omega + 0.1) + (x1_dot(3) - uCR*vL)^2 - z(2)]; % isolated nonlinearity
n_out = length(phi);
% concatination of all states&controls Phi is dependent on
y = vertcat(xL, vL, theta, omega, x1_dot(3), z(2));
uhat = uCR;
ny = length(y);
nuhat = length(uhat);

C = zeros( nx1+nz, n_out); C(8,1) = 1; C(9,2) = 1; C(10,3) = 1;

% linear input matrices
L_x_fun = Function('L_x_fun',{x1},{jacobian(y,x1)});
L_xdot_fun = Function('L_x_fun',{x1},{jacobian(y,x1_dot)});
L_z_fun = Function('L_z_fun',{x1},{jacobian(y,z)});

L_u_fun = Function('L_u_fun',{x1},{jacobian(uhat,u)});

L_x = full(L_x_fun(0));
L_xdot = full(L_xdot_fun(0));
L_u = full(L_u_fun(0));
L_z = full(L_z_fun(0));

y_check = L_xdot * x1_dot +L_x * x1 + L_z * z; %% THis should be the same as y
uhat_check = L_u * u;

jac_phi_y = jacobian(phi,y);
jac_phi_uhat = jacobian(phi,uhat);

phi_fun = Function([model_name_prefix,'phi_fun'], {y,uhat}, {phi});
phi_fun_jac_y = Function([model_name_prefix,'phi_fun_jac_y'], {y,uhat}, {phi, jac_phi_y});
phi_jac_y_uhat = Function([model_name_prefix,'phi_jac_y_uhat'], {y,uhat}, {jac_phi_y, jac_phi_uhat});

phi_jac_y = Function([model_name_prefix,'phi_jac_y_uhat'], {y,uhat}, {jac_phi_y});

% Linear output
ALO = zeros(nx2);
% A2(1,1) = 1;

% f = uCR^2 + xL^2;
jac_f_x1 = jacobian(f,x1);
jac_f_u = jacobian(f,u);
jac_f_z = jacobian(f,z);
jac_f_k1 = jacobian(f,x1_dot);

f_fun = Function('f_los', {x1_dot,x1,z,u}, {f});

% jac_Phi_u_fun = Function('jac_Phi_u_fun', {y,u},{jac_Phi_u});

f_lo_fun_jac_x1k1uz = Function([model_name_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1_dot, z, u}, ...
    {f, [jac_f_x1, jac_f_k1, jac_f_z, jac_f_u]});

f_lo_fun = Function([model_name_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1_dot, z, u}, {f});
% struct for matlab prototype
s = struct('A', A, 'B', B, 'C', C, 'E', E, 'ALO',ALO, 'L_x', L_x, 'L_xdot', L_xdot, 'L_z', L_z, 'L_u', L_u, ...
    'phi_fun_jac_y', phi_fun_jac_y, 'phi_jac_y_uhat', phi_jac_y_uhat, 'f_fun', f_fun, ...
    'nx1', nx1, 'nx2', nx2, 'nu', nu, 'n_out', n_out, 'nx', nx, 'nz', nz, 'ny', ny, 'nuhat', nuhat,...
    'f_lo_fun_jac_x1k1uz', f_lo_fun_jac_x1k1uz);


%% generate functions
% get matrices
dummy = SX.sym('dummy');

model_matrices = SX.zeros(size([A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:); c(:)])) + ...
    [A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:); c(:)];
get_matrices_fun = Function([model_name_prefix,'get_matrices_fun'], {dummy}, {model_matrices(:)});
get_matrices_fun.generate([model_name_prefix,'get_matrices_fun'], casadi_opts);

% generate Phi, f_LO
f_lo_fun_jac_x1k1uz.generate([model_name_prefix,'f_lo_fun_jac_x1k1uz'], casadi_opts);
phi_fun.generate([model_name_prefix,'phi_fun'], casadi_opts);
phi_fun_jac_y.generate([model_name_prefix,'phi_fun_jac_y'], casadi_opts);
phi_jac_y_uhat.generate([model_name_prefix,'phi_jac_y_uhat'], casadi_opts);

% %% generate functions
% casadi_opts = struct('mex', false);
% dummy = SX.sym('dummy');
% 
% % get ints
n_steps = 1;
num_stages = 4;

%% test model function
x0 = zeros(nx,1);
x0(3) = .8;
x0dot = zeros(nx,1);
z0 = zeros(nz,1);
u0 = [40.108149413030752; 50.446662212534974];

x1_0 = x0(1:nx1);
x1dot_0 = x0dot(1:nx1);
x2_0 = x0(nx1+1:nx);

y0 = L_x * x1_0 + L_xdot * x1dot_0 + L_z * z0;

uhat0 = L_u * u0;
phi_val0 = phi_fun(y0, uhat0);
gnsf_val0 = [E * [x1dot_0; z0] - full(A*x1_0 + B * u0 + C * phi_val0 + c); 
            full(f_lo_fun(x1_0, x1dot_0, z0, u0))]

%% set up equivalent implicit model
x2_dot = SX.sym('x2_dot', nx2);
xdot = [x1_dot; x2_dot];
f_impl = [(f_expl - [x1_dot; x2_dot]);
            z - [((theta^2)/8 + xL);
      cos(omega + 0.1) + (x1_dot(3) - uCR*vL)^2 - z(2)]];

% impl_ode_fun
impl_ode_fun = Function([model_name_prefix,'impl_ode_fun'], {x, xdot, u, z}, {f_impl});
impl_ode_fun.generate([model_name_prefix,'impl_ode_fun'], casadi_opts);

% impl_ode_fun_jac_x_xdot
impl_ode_fun_jac_x_xdot = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot'], {x, xdot, u, z}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});
impl_ode_fun_jac_x_xdot.generate([model_name_prefix,'impl_ode_fun_jac_x_xdot'], casadi_opts);

% impl_ode_jac_x_xdot_u
impl_ode_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_jac_x_xdot_u'], {x, xdot, u, z}, {jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, u), jacobian(f_impl, z)});
impl_ode_jac_x_xdot_u.generate([model_name_prefix,'impl_ode_jac_x_xdot_u'], casadi_opts);

% impl_fun_ode_jac_x_xdot_u
impl_ode_fun_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot_u'], {x, xdot, u, z}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});
impl_ode_fun_jac_x_xdot_u.generate([model_name_prefix,'impl_ode_fun_jac_x_xdot_u'], casadi_opts);

%% check if equivalent
f_imp_0 = full(impl_ode_fun(x0, x0dot, u0, z0))

        




 
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
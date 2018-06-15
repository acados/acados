clc;
clear;
close all;

%% CasADi
import casadi.*
if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	error('Please download and install Casadi 3.4.0')
end
model_name_prefix = 'inv_pendulum_';

%% Parameters (taken from Rien Quirynens Master Thesis)
m = 2;
g = 9.81;
M = 3.5;
I = 0.1;

%% Set up States & Controls
xpos    = SX.sym('xpos');     % Differential States
ypos    = SX.sym('ypos');
alpha   = SX.sym('alpha');     
vx      = SX.sym('vx');
vy      = SX.sym('vy');
valpha  = SX.sym('valpha');
x = vertcat(xpos, ypos, alpha, vx, vy, valpha);

ax      = SX.sym('ax');     % Algebraic states
ay      = SX.sym('ay');
aalpha  = SX.sym('aalpha');
Fx      = SX.sym('Fx');
Fy      = SX.sym('Fy');
z = vertcat(ax, ay, aalpha, Fx, Fy);

u       = SX.sym('u');  % Controls

%% xdot
xpos_dot    = SX.sym('xpos_dot');     % Differential States
ypos_dot    = SX.sym('ypos_dot');
alpha_dot   = SX.sym('alpha_dot');     
vx_dot      = SX.sym('vx_dot');
vy_dot      = SX.sym('vy_dot');
valpha_dot  = SX.sym('valpha_dot');

xdot = [xpos_dot; ypos_dot; alpha_dot; vx_dot; vy_dot; valpha_dot];

%% dimensions
nx = length(x);
nu = length(u);
nz = length(z);

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

%% initial value
x0 = [1; -5; 1; 0.1; -0.5; 0.1];
z0 = [-1.5; -0.3; -0.3; -3; 19];
u0 = 1;

%% Generalized nonlinear static feedback formulation (GNSF)

x1 = x;
x1_dot = xdot;
nx1 = length(x1);
x2 = []; % Linear output/quadrature state
f = [];
nx2 = 0;

%% Model defining matrices
E = zeros(nx1 + nz);
E(1:nx1, 1:nx1) = eye(nx1); % equs corresponding to x1 are explicit

A = zeros(nx1 + nz, nx1);
B = zeros(nx1 + nz, nu);
c = zeros(nx1 + nz, 1);
phi = [];
% write linear terms into corresponding matrices
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;

E(4, nx1 + 1) = -1;
E(5, nx1 + 2) = -1;
E(6, nx1 + 3) = -1;

E(7, nx1 + 1) = m;
E(7, nx1 + 4) = -1;
B(7,1) = 1;

% m * ay + m * g - Fy
E(8, nx1 + 2) = m;
E(8, nx1 + 5) = -1;
c(8) = -m * g;

% I * aalpha - M - (Fx+u)*y + Fy * x
E(9, nx1 + 3) = I;
c(9) = M;
phi = [phi; (Fx + u) * ypos - Fy * xpos];
C(9,1) = 1;

% ax + vy * valpha + ypos * aalpha
E(10, nx1 + 1) = 1;
phi = [phi; -(vy * valpha + ypos * aalpha)];
C(10,2) = 1;

% ay - vx * valpha - xpos * aalpha;
E(11, nx1 + 2) = 1;
phi = [phi; (vx * valpha + xpos * aalpha)];
C(11, 3) = 1;

% check invertibility
E11 = E(1:nx1, 1:nx1);
E22 = E(nx1+1:nx1+nz, nx1+1:nx1+nz);

assert(rank(E11) == nx1)
assert(rank(E22) == nz)

%% regard phi and collect differential & algebraic states it depends on in y
y = [xpos; ypos; vx; vy; valpha; aalpha; Fx; Fy];
% egard phi and collect controls it depends on in uhat
uhat = u;

ny = length(y);
nuhat = length(uhat);
n_out = length(phi);


% linear input matrices
L_x_fun = Function('L_x_fun',{x1},{jacobian(y,x1)});
L_xdot_fun = Function('L_x_fun',{x1},{jacobian(y,x1_dot)});
L_z_fun = Function('L_z_fun',{x1},{jacobian(y,z)});
L_u_fun = Function('L_u_fun',{x1},{jacobian(uhat,u)});

L_x = full(L_x_fun(0));
L_xdot = full(L_xdot_fun(0));
L_u = full(L_u_fun(0));
L_z = full(L_z_fun(0));

%% generate nonlinearity function phi
jac_phi_y = jacobian(phi,y);
jac_phi_uhat = jacobian(phi,uhat);

phi_fun = Function([model_name_prefix,'phi_fun'], {y,uhat}, {phi});
phi_fun_jac_y = Function([model_name_prefix,'phi_fun_jac_y'], {y,uhat}, {phi, jac_phi_y});
phi_jac_y_uhat = Function([model_name_prefix,'phi_jac_y_uhat'], {y,uhat}, {jac_phi_y, jac_phi_uhat});
phi_jac_y = Function([model_name_prefix,'phi_jac_y_uhat'], {y,uhat}, {jac_phi_y});


%% Linear output system
nx2 = 0;
ALO = [];
f = [];

jac_f_x1 = jacobian(f,x1);
jac_f_u = jacobian(f,u);
jac_f_z = jacobian(f,z);
jac_f_k1 = jacobian(f,x1_dot);

% generate
f_lo_fun_jac_x1k1uz = Function([model_name_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1_dot, z, u}, ...
    {f, [jac_f_x1, jac_f_k1, jac_f_z, jac_f_u]});
f_lo_fun = Function([model_name_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1_dot, z, u}, ...
    {f}); % just to check model

%% generate get matrices
dummy = SX.sym('dummy');
model_matrices = [A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:); c(:)];
get_matrices_fun = Function([model_name_prefix,'get_matrices_fun'], {dummy}, {model_matrices(:)});

%% export implicit model
% impl_ode_fun
impl_ode_fun = Function([model_name_prefix,'impl_ode_fun'],...
                        {x, xdot, u, z}, {f_impl});

% impl_ode_fun_jac_x_xdot
impl_ode_fun_jac_x_xdot = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot'],...
          {x, xdot, u, z}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});

% impl_ode_jac_x_xdot_u
impl_ode_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_jac_x_xdot_u'],...
    {x, xdot, u, z}, {jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, u), jacobian(f_impl, z)});

% impl_fun_ode_jac_x_xdot_u
impl_ode_fun_jac_x_xdot_u = Function([model_name_prefix,'impl_ode_fun_jac_x_xdot_u'],...
    {x, xdot, u, z}, {f_impl, jacobian(f_impl, x), jacobian(f_impl, xdot), jacobian(f_impl, z)});

%% test model function -- check equivalence
% x0dot = zeros(nx,1);
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

%% Generate C functions
EXPORT = 1;

if EXPORT
    impl_ode_fun.generate([model_name_prefix,'impl_ode_fun'], casadi_opts);
    impl_ode_fun_jac_x_xdot.generate([model_name_prefix,'impl_ode_fun_jac_x_xdot'], casadi_opts);
    impl_ode_jac_x_xdot_u.generate([model_name_prefix,'impl_ode_jac_x_xdot_u'], casadi_opts);
    impl_ode_fun_jac_x_xdot_u.generate([model_name_prefix,'impl_ode_fun_jac_x_xdot_u'], casadi_opts);
    phi_fun.generate([model_name_prefix,'phi_fun'], casadi_opts);
    phi_fun_jac_y.generate([model_name_prefix,'phi_fun_jac_y'], casadi_opts);
    phi_jac_y_uhat.generate([model_name_prefix,'phi_jac_y_uhat'], casadi_opts);
    get_matrices_fun.generate([model_name_prefix,'get_matrices_fun'], casadi_opts);
    f_lo_fun_jac_x1k1uz.generate([model_name_prefix,'f_lo_fun_jac_x1k1uz'], casadi_opts);
end

% test jac_x
% phi_jac_y0 = full(phi_jac_y(y0, uhat0)) * L_x

       




 
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
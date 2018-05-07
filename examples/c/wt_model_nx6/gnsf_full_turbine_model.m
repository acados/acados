clear all;
close all;
clc
%restoredefaultpath

%% perpare the MATLAB path to include CASADI
%tmp = path();
%if(isempty(strfind(strrep(tmp,'\','/'),'D:\temp_Axel\casadi-mat2014')))
%    addpath(path,'D:\temp_Axel\casadi-mat2014')
%end;

import casadi.*

%% define the symbolic variables of the plant
S02_DefACADOSVarSpace_gnsf;

%% load plant parameters
S03_SetupSysParameters;

%% Define casadi spline functions
% aerodynamic torque coefficient for FAST 5MW reference turbine
load('CmDataSpline.mat')
c_StVek = c_St';
splineCMBL = interpolant('Spline','bspline',{y_St,x_St},c_StVek(:));
clear x_St y_St c_St c_StVek

%% define ode rhs in explicit form (22 equations)
S04_SetupNonlinearStateSpaceDynamics;


%% Generalized nonlinear static feedback formulation (GNSF)
casadi_opts = struct('mex', false);
x = [xy; xnoy];

nx = length(x);
nu = length(u);
np = length(p);

x = [xy; xnoy];
x1 = x;
nx1 = length(x1);
z = MX.sym('z',0);
nz = 0;
% x2 = SX.sy('x2',0);
nx2 = 0;
x1_dot = MX.sym('x1_dot',nx1,1);

% if CasadiMeta.version()=='3.4.0'
% 	% casadi 3.4
% 	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
% else
% 	% old casadi versions
% 	casadi_opts = struct('mex', false);
% end
casadi_export_prefix = 'casadi_';
casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
casadi_opts_mex = struct('mex', true, 'casadi_int', 'int', 'casadi_real', 'double');

%% Model defining matrices
A = zeros(nx);
A(1,4) = p_14/(p_10+p_11);
A(1,2) = p_13/(p_10+p_11);
A(1,6) = -p_12/(p_10+p_11);
A(3,1) = -p_8;
A(4,2) = -p_8;
A(5,5) = -p_15;

B = zeros(nx, nu);
B(5,1) = p_15;
A(6,6) = -p_16;
B(6,2) = p_16;

phi = fe(2);

n_out  = length(phi);
C = zeros(nx, n_out); C(2,1) = 1;

E = eye(nx+nz);
y = [x(1:5)];


uhat = [];
ny = length(y);
nuhat = length(uhat);

% linear input matrices
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

phi_fun = Function([casadi_export_prefix,'phi_fun'], {y,uhat,p}, {phi});
phi_fun_jac_y = Function([casadi_export_prefix,'phi_fun_jac_y'], {y,uhat,p}, {phi, jac_phi_y});
phi_jac_y_uhat = Function([casadi_export_prefix,'phi_jac_y_uhat'], {y,uhat,p}, {jac_phi_y, jac_phi_uhat});

phi_jac_y = Function([casadi_export_prefix,'phi_jac_y_uhat'], {y,uhat,p}, {[jac_phi_y]});

% Linear output
ALO = zeros(nx2);
% A2(1,1) = 1;

f = [];
% f = uCR^2 + xL^2;
jac_f_x1 = jacobian(f,x1);
jac_f_u = jacobian(f,u);
jac_f_z = jacobian(f,z);
jac_f_k1 = jacobian(f,x1_dot);

f_fun = Function('f_los', {x1_dot,x1,z,u}, {f});

% jac_Phi_u_fun = Function('jac_Phi_u_fun', {y,u},{jac_Phi_u});

f_lo_fun_jac_x1k1uz = Function([casadi_export_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1_dot, z, u}, ...
    {f, [jac_f_x1, jac_f_k1, jac_f_z, jac_f_u]});

% struct for matlab prototype
s = struct('A', A, 'B', B, 'C', C, 'E', E, 'ALO',ALO, 'L_x', L_x, 'L_xdot', L_xdot, 'L_z', L_z, 'L_u', L_u, ...
    'phi_fun_jac_y', phi_fun_jac_y, 'phi_jac_y_uhat', phi_jac_y_uhat, 'f_fun', f_fun, ...
    'nx1', nx1, 'nx2', nx2, 'nu', nu, 'n_out', n_out, 'nx', nx, 'nz', nz, 'ny', ny, 'nuhat', nuhat,...
    'f_lo_fun_jac_x1k1uz', f_lo_fun_jac_x1k1uz);


%% generate functions
% ints = SX.zeros(8,1) + [s.nx, s.nu, s.nz, s.nx1, s.nx2, q, n_steps, s.n_out]';
% get_ints_fun = Function('get_ints_fun',{x},{[s.nx, s.nu, s.nz, s.nx1, s.nx2, q, n_steps, s.n_out]});
%     get_ints_fun.generate('get_ints_fun', casadi_opts);

% get matrices
dummy = SX.sym('dummy');

model_matrices = SX.zeros(size([A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:)])) + ...
    [A(:); B(:); C(:); E(:); L_x(:); L_xdot(:); L_z(:); L_u(:); ALO(:)];
get_matrices_fun = Function([casadi_export_prefix,'get_matrices_fun'], {dummy}, {model_matrices(:)});
get_matrices_fun.generate('get_matrices_fun', casadi_opts);

% generate Phi, f_LO
f_lo_fun_jac_x1k1uz.generate(['f_lo_fun_jac_x1k1uz'], casadi_opts);
phi_fun.generate(['phi_fun'], casadi_opts);
phi_fun_jac_y.generate(['phi_fun_jac_y'], casadi_opts);
phi_jac_y_uhat.generate(['phi_jac_y_uhat'], casadi_opts);
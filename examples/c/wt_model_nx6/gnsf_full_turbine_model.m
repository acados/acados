%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

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

if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
    error('Provide Casadi version 3.4.0');
end
casadi_export_prefix = 'casadi_';

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

c = zeros(nx1 + nz,1);

phi = fe(2);

n_out  = length(phi);
C = zeros(nx, n_out); C(2,1) = 1;

E = eye(nx+nz);
y = [x(1:5)];


uhat = [];
ny = length(y);
nuhat = length(uhat);

% linear input matrices
L_x_fun = Function('L_x_fun',{x1},{jacobian(y,x1)});
L_xdot_fun = Function('L_x_fun',{x1},{jacobian(y,x1_dot)});
L_z_fun = Function('L_z_fun',{x1},{jacobian(y,z)});

L_u_fun = Function('L_u_fun',{x1},{jacobian(uhat,u)});

L_x = full(L_x_fun(0));
L_xdot = full(L_xdot_fun(0));
L_u = full(L_u_fun(0));
L_z = full(L_z_fun(0));

y_check = L_xdot * x1_dot +L_x * x1 + L_z * z; % This should be the same as y
uhat_check = L_u * u;

jac_phi_y = jacobian(phi,y);
jac_phi_uhat = jacobian(phi,uhat);

phi_fun = Function([casadi_export_prefix,'phi_fun'], {y,uhat,p}, {phi});
phi_fun_jac_y = Function([casadi_export_prefix,'phi_fun_jac_y'], {y,uhat,p}, {phi, jac_phi_y});
phi_jac_y_uhat = Function([casadi_export_prefix,'phi_jac_y_uhat'], {y,uhat,p}, {jac_phi_y, jac_phi_uhat});

phi_jac_y = Function([casadi_export_prefix,'phi_jac_y_uhat'], {y,uhat,p}, {[jac_phi_y]});

% Linear output
A_LO = zeros(nx2);
E_LO = eye(nx2);
% A2(1,1) = 1;

f = [];
% f = uCR^2 + xL^2;
jac_f_x1 = jacobian(f,x1);
jac_f_u = jacobian(f,u);
jac_f_z = jacobian(f,z);
jac_f_k1 = jacobian(f,x1_dot);

f_fun = Function('f_los', {x1_dot,x1,z,u}, {f});

f_lo_fun_jac_x1k1uz = Function([casadi_export_prefix,'f_lo_fun_jac_x1k1uz'], {x1, x1_dot, z, u}, ...
    {f, [jac_f_x1, jac_f_k1, jac_f_z, jac_f_u]});


%% generate functions
% get matrices
dummy = SX.sym('dummy');

get_matrices_fun = Function([casadi_export_prefix,'get_matrices_fun'], {dummy},...
     {A, B, C, E, L_x, L_xdot, L_z, L_u, A_LO, c, E_LO});
get_matrices_fun.generate('get_matrices_fun', casadi_opts);

% generate Phi, f_LO
f_lo_fun_jac_x1k1uz.generate(['f_lo_fun_jac_x1k1uz'], casadi_opts);
phi_fun.generate(['phi_fun'], casadi_opts);
phi_fun_jac_y.generate(['phi_fun_jac_y'], casadi_opts);
phi_jac_y_uhat.generate(['phi_jac_y_uhat'], casadi_opts);


%% to get value for IRK initialization equivalent to the GNSF 0 init.
x0 = [1.353969828015453e+00, 3.228986350219792e-04, 1.799917898419573e+02, 4.306609995588871e-03, 8.269214720965625e+00, 4.010510098723469e+00]';
% u0 = [8.169651470932033e+00, 4.024634365037572e+00]';
u0 = [0;0];
p0 = 1;

format long e
disp('value for equivalent irk initialization')
xdot1_z1_0 = (E\(A*x0(1:nx1)+ B * u0 + c))'


%% check if same result as in full_turbine_model.m
disp('value to check if result is the same as in full_turbine_model.m')
x0 = ones(nx,1);
x0dot = ones(nx,1);
u0 = ones(nu,1);
p0 = 1;

y0 = L_x * x0 + L_xdot * x0dot;
uhat0 = L_u * u0;

gnsf0 = E*x0dot - A*x0 - B*u0 - C*phi_fun(y0, uhat0, p0) - c

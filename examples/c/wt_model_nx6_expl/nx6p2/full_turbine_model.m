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
S02_DefACADOSVarSpace;

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

%% generate casadi C functions
nx = 8;
nu = 3;

% ODE
odeFun = Function('expl_ode', {x, u}, {f});
odeFun.generate('expl_ode');

% jac x
jac_x_odeFun = Function('expl_ode_jac_x', {x, u}, {jacobian(f, x)});
jac_x_odeFun.generate('expl_ode_jac_x');

% jac u
jac_u_odeFun = Function('expl_ode_jac_u', {x, u}, {jacobian(f, u)});
jac_u_odeFun.generate('expl_ode_jac_u');

% forward VDE
Sx = MX.sym('Sx', nx, nx);
Su = MX.sym('Su', nx, nu);

%vdeX = MX.zeros(nx, nx) + jtimes(f, x, Sx);
vdeX = MX.zeros(nx, nx) + jacobian(f, x)*Sx;

%vdeU = MX.zeros(nx, nu) + jtimes(f, x, Su) + jacobian(f, u);
vdeU = MX.zeros(nx, nu) + jacobian(f, x)*Su + jacobian(f, u);

forw_vdeFun = Function('expl_forw_vde', {x, Sx, Su, u}, {f, vdeX, vdeU});
forw_vdeFun.generate('expl_forw_vde');


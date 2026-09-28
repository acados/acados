%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function model = ocp_model_wind_turbine_nx6()

import casadi.*

%% define the symbolic variables of the plant
ocp_S02_DefACADOSVarSpace;

%% load plant parameters
ocp_S03_SetupSysParameters;

%% Define casadi spline functions
% aerodynamic torque coefficient for FAST 5MW reference turbine
load('CmDataSpline.mat')
c_StVek = c_St';

% set different bspline degree default is cubic
% opt = struct;
% opt = [1; 1];
% opt = [3; 3];
% opt = [5; 5];
% splineCMBL = interpolant('Spline','bspline',{y_St,x_St},c_StVek(:), opt);

splineCMBL = interpolant('Spline','bspline',{y_St,x_St},c_StVek(:));
clear x_St y_St c_St c_StVek

%% define ode rhs in explicit form (22 equations)
ocp_S04_SetupNonlinearStateSpaceDynamics;

%% generate casadi C functions
nx = 8;
nu = 2;
np = 1;

%% populate structure
model = AcadosModel();
model.name = 'wind_turbine_nx6';
model.x = x;
model.u = u;
model.xdot = dx;
model.p = p;
model.f_expl_expr = f_expl;
model.f_impl_expr = f_impl;
model.con_h_expr = h;
model.con_h_expr_e = h_e;

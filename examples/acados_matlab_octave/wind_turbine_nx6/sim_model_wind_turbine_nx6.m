%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



% NOTE: `acados` currently supports both an old MATLAB/Octave interface (< v0.4.0)
% as well as a new interface (>= v0.4.0).

% THIS EXAMPLE still uses the OLD interface. If you are new to `acados` please start
% with the examples that have been ported to the new interface already.
% see https://github.com/acados/acados/issues/1196#issuecomment-2311822122)


function model = sim_model_wind_turbine_nx6()

import casadi.*

%% define the symbolic variables of the plant
sim_S02_DefACADOSVarSpace;

%% load plant parameters
sim_S03_SetupSysParameters;

%% Define casadi spline functions
% aerodynamic torque coefficient for FAST 5MW reference turbine
load('CmDataSpline.mat')
c_StVek = c_St';

% bspline by CasADi default are cubic
splineCMBL = interpolant('Spline','bspline',{y_St,x_St},c_StVek(:));
clear x_St y_St c_St c_StVek

%% define ode rhs in explicit form (22 equations)
sim_S04_SetupNonlinearStateSpaceDynamics;

%% generate casadi C functions
nx = 6;
nu = 2;
np = 1;

%% populate structure
model = AcadosModel();
model.name = 'wind_turbine_nx6';
model.x = x;
model.u = u;
model.xdot = dx;
model.p = p;
model.f_expl_expr = fe;
model.f_impl_expr = f_impl;
%model.expr_h = h;
%model.expr_h_e = hN;
%model.expr_y = expr_y;
%model.expr_y_e = expr_y_e;

end
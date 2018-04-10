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


% casadi opts for code generation
if CasadiMeta.version()=='3.4.0'
	% casadi 3.4
	opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
else
	% old casadi versions
	opts = struct('mex', false);
end


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
nu = 2;

% expl_ode_fun

expl_ode_fun = Function('casadi_expl_ode_fun', {x, u, p}, {fe});
expl_ode_fun.generate('expl_ode_fun', opts);


% expl_vde_for

Sx = MX.sym('Sx', nx, nx);
Su = MX.sym('Su', nx, nu);

vdeX = MX.zeros(nx, nx) + jtimes(fe, x, Sx);
%vdeX = MX.zeros(nx, nx) + jacobian(fe, x)*Sx;

vdeU = MX.zeros(nx, nu) + jtimes(fe, x, Su) + jacobian(fe, u);
%vdeU = MX.zeros(nx, nu) + jacobian(fe, x)*Su + jacobian(fe, u);

expl_vde_for = Function('casadi_expl_vde_for', {x, Sx, Su, u, p}, {fe, vdeX, vdeU});
expl_vde_for.generate('expl_vde_for', opts);


% impl_ode_fun

impl_ode_fun = Function('casadi_impl_ode_fun', {x, dx, u, p}, {fi});
impl_ode_fun.generate('impl_ode_fun', opts);


% impl_ode_jac_x

impl_ode_jac_x = Function('casadi_impl_ode_jac_x', {x, dx, u, p}, {jacobian(fi, x)});
impl_ode_jac_x.generate('impl_ode_jac_x', opts);


% impl_ode_jac_xdot

impl_ode_jac_xdot = Function('casadi_impl_ode_jac_xdot', {x, dx, u, p}, {jacobian(fi, dx)});
impl_ode_jac_xdot.generate('impl_ode_jac_xdot', opts);


% impl_ode_jac_u

impl_ode_jac_u = Function('casadi_impl_ode_jac_u', {x, dx, u, p}, {jacobian(fi, u)});
impl_ode_jac_u.generate('impl_ode_jac_u', opts);


% impl_ode_fun_jac_x_xdot

impl_ode_fun_jac_x_xdot = Function('casadi_impl_ode_fun_jac_x_xdot', {x, dx, u, p}, {fi, jacobian(fi, x), jacobian(fi, dx)});
impl_ode_fun_jac_x_xdot.generate('impl_ode_fun_jac_x_xdot', opts);


% impl_ode_jac_x_xdot_u

impl_ode_jac_x_xdot_u = Function('casadi_impl_ode_jac_x_xdot_u', {x, dx, u, p}, {jacobian(fi, x), jacobian(fi, dx), jacobian(fi, u)});
impl_ode_jac_x_xdot_u.generate('impl_ode_jac_x_xdot_u', opts);


% impl_ode_jac_x_u

impl_ode_jac_x_u = Function('casadi_impl_ode_jac_x_u', {x, dx, u, p}, {jacobian(fi, x), jacobian(fi, u)});
impl_ode_jac_x_u.generate('impl_ode_jac_x_u', opts);




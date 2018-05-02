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
	error('Please download and install Casadi 3.4.0')
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
nx = 6;
nu = 2;


% expl_ode_fun

expl_ode_fun = Function('casadi_expl_ode_fun', {x, u, p}, {fe});
expl_ode_fun.generate('expl_ode_fun', opts);


% expl_vde_for

Sx = MX.sym('Sx', nx, nx);
Su = MX.sym('Su', nx, nu);

%vdeX = MX.zeros(nx, nx) + jtimes(fe, x, Sx);
vdeX = MX.zeros(nx, nx) + jacobian(fe, x)*Sx;

%vdeU = MX.zeros(nx, nu) + jtimes(fe, x, Su) + jacobian(fe, u);
vdeU = MX.zeros(nx, nu) + jacobian(fe, x)*Su + jacobian(fe, u);

expl_vde_for = Function('casadi_expl_vde_for', {x, Sx, Su, u, p}, {fe, vdeX, vdeU});
expl_vde_for.generate('expl_vde_for', opts);


% expl_vde_adj

lam = MX.sym('lam', nx, 1);

adj = jtimes(fe, [x; u], lam, true);

expl_vde_adj = Function('casadi_expl_vde_adj', {x, lam, u, p}, {adj});
expl_vde_adj.generate('expl_vde_adj', opts);


% impl_ode_fun

impl_ode_fun = Function('casadi_impl_ode_fun', {x, dx, u, p}, {fi});
impl_ode_fun.generate('impl_ode_fun', opts);

% impl_ode_fun_jac_x_xdot

impl_ode_fun_jac_x_xdot = Function('casadi_impl_ode_fun_jac_x_xdot', {x, dx, u, p}, {fi, jacobian(fi, x), jacobian(fi, dx)});
impl_ode_fun_jac_x_xdot.generate('impl_ode_fun_jac_x_xdot', opts);


% impl_ode_jac_x_xdot_u

impl_ode_jac_x_xdot_u = Function('casadi_impl_ode_jac_x_xdot_u', {x, dx, u, p}, {jacobian(fi, x), jacobian(fi, dx), jacobian(fi, u)});
impl_ode_jac_x_xdot_u.generate('impl_ode_jac_x_xdot_u', opts);


% impl_ode_fun_jac_x_xdot_u

impl_ode_fun_jac_x_xdot_u = Function('casadi_impl_ode_fun_jac_x_xdot_u', {x, dx, u, p}, {fi, jacobian(fi, x), jacobian(fi, dx), jacobian(fi, u)});
impl_ode_fun_jac_x_xdot_u.generate('impl_ode_fun_jac_x_xdot_u', opts);

return







%% create an ODE casadi object
ode = struct('x',x,'p',u,'ode',fe);
%% Instantiate casadi integrator object -> using fixed-step Runge Kutta of order 4
Ts = 0.2;
nbrIntermedSamples = 10;
ts = linspace(0,Ts,nbrIntermedSamples);         % time grid for each integration step
opts = struct('grid', ts,'output_t0', 1,'print_stats',1);
casadiIntObj = casadi.integrator('I', 'rk', ode, opts);
%% load reference data for simulation
load('testSim.mat')
x0       = [statesFAST(1,:)];     % initial state for starting the simulation

%% simulate dynamics in a step-wise fashion
len = length(tFAST);
xTraj = x0;     % storage element for simulated state trajectory

% to avoid unstable behavior introduce a small pi-controller for rotor
% speed tracking
uctrl = 0;
uctrlI = 0;
kI = 1e-1;
kP = 10;
Ck = [];
for ii=1:len-1
    
    % compile inputs (parameters) for current step
    u0 = [Usim(ii,:)];
    u0(2) = max(u0(2) - uctrl,0);
    
    % display simulation progess
    if(mod(ii,10)==0)
        display(['Simulation time t = ' num2str(tFAST(ii)) ' ...']);
    end;
    
    % execute simulation step with current input and state
    res = casadiIntObj('x0', x0, 'p', u0);
                    
    % extract state at next time step
    xTraj = cat(1,xTraj,full(res.xf(:,end))');
    % update initial state for subsequent simulation step
    x0 = xTraj(end,:);
    
    % update PI-controller
    ctrlErr = statesFAST(ii+1,1)-x0(1);
    uctrlI = uctrlI + kI*ctrlErr*Ts;
    uctrl = kP*ctrlErr + uctrlI;
end;

%% Plot the simulation results and compare to reference simulation
close all;
x_output = xTraj';
% plot all turbine states
for ii=1:4:4
    figure
    subplot(4,1,1)
    plot(tFAST(1:len-1),x_output(ii,1:len-1),tFAST(1:len-1),statesFAST(1:len-1,ii))
    legend({'casadi','FAST'})
    subplot(4,1,2)
    plot(tFAST(1:len-1),x_output(ii+1,1:len-1),tFAST(1:len-1),statesFAST(1:len-1,ii+1))
    legend({'casadi','FAST'})
    subplot(4,1,3)
    plot(tFAST(1:len-1),x_output(ii+2,1:len-1),tFAST(1:len-1),statesFAST(1:len-1,ii+2))
    legend({'casadi','FAST'})
    subplot(4,1,4)
    plot(tFAST(1:len-1),x_output(ii+3,1:len-1),tFAST(1:len-1),statesFAST(1:len-1,ii+3))
    legend({'casadi','FAST'})
end;

% plot actuator states
figure
subplot(2,1,1)
plot(tFAST(1:len-1),x_output(5,1:len-1))
subplot(2,1,2)
plot(tFAST(1:len-1),x_output(6,1:len-1))

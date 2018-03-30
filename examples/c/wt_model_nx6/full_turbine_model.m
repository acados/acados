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
nx = 6;
nu = 2;

% ODE
odeFun = Function('casadi_expl_ode_fun', {x, u, p}, {f});
odeFun.generate('expl_ode_fun');

% jac x
jac_x_odeFun = Function('casadi_expl_ode_jac_x', {x, u, p}, {jacobian(f, x)});
jac_x_odeFun.generate('expl_ode_jac_x');

% jac u
jac_u_odeFun = Function('casadi_expl_ode_jac_u', {x, u, p}, {jacobian(f, u)});
jac_u_odeFun.generate('expl_ode_jac_u');

% forward VDE
Sx = MX.sym('Sx', nx, nx);
Su = MX.sym('Su', nx, nu);

%vdeX = MX.zeros(nx, nx) + jtimes(f, x, Sx);
vdeX = MX.zeros(nx, nx) + jacobian(f, x)*Sx;

%vdeU = MX.zeros(nx, nu) + jtimes(f, x, Su) + jacobian(f, u);
vdeU = MX.zeros(nx, nu) + jacobian(f, x)*Su + jacobian(f, u);

forw_vdeFun = Function('casadi_expl_vde_for', {x, Sx, Su, u, p}, {f, vdeX, vdeU});
forw_vdeFun.generate('expl_vde_for');

return







%% create an ODE casadi object
ode = struct('x',x,'p',u,'ode',f);
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

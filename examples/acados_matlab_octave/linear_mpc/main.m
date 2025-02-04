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

clear all; clc; close all;
check_acados_requirements()

%% solver settings
h = 1/14;               % [s] sampling time
N_horizon = 10;         % [-] number of prediction steps
tf = N_horizon * h;     % [s] prediction horizon length

[model,Ad,Bd] = quadcopter_model();     % system dynamics
nx = length(model.x);                   % state size
nu = length(model.u);                   % input size

% ocp model
ocp = AcadosOcp();
ocp.name = 'quad_mpc';
ocp.model = model;

%% linear least squares cost
% cost matrices
Q = diag([0 0 10 10 10 10 0 0 0 5 5 5]);    % state cost
R = 0.1*eye(nu);                            % input cost

% details on discrete cost scaling:
% https://docs.acados.org/python_interface/index.html#acados_template.acados_ocp_cost.AcadosOcpCost

% initial cost (inputs only)
ocp.cost.cost_type_0 = 'LINEAR_LS';
ny_0 = nu;
ocp.cost.Vu_0 = eye(ny_0);
ocp.cost.Vx_0 = zeros(ny_0,nx);
ocp.cost.W_0 = 1/h * 2 * R;  % scale to match the original cost
ocp.cost.yref_0 = zeros(ny_0,1);

% path cost (inputs and states)
ocp.cost.cost_type = 'LINEAR_LS';
ny = nu + nx;
ocp.cost.Vu = [eye(nu); zeros(nx,nu)];
ocp.cost.Vx = [zeros(nu, nx); eye(nx)];
ocp.cost.W = 1/h * 2 * blkdiag(R,Q);  % scale to match the original cost
ocp.cost.yref = zeros(ny, 1);

% terminal cost (states only)
ocp.cost.cost_type_e = 'LINEAR_LS';
ny_e = nx;
ocp.cost.Vx_e = eye(ny_e, nx);
ocp.cost.W_e = 2 * Q;  % terminal cost is not scaled with h later
ocp.cost.yref_e = zeros(ny_e, 1);

%% constraints
% input constraints
u0 = 10.5916;                                       % steady-state input
ocp.constraints.idxbu = 0:3;                        % zero-based indices
ocp.constraints.lbu = [9.6; 9.6; 9.6; 9.6] - u0;    % input lower bounds
ocp.constraints.ubu = [13; 13; 13; 13] - u0;        % input upper bounds

% placeholder for the initial state constraint
ocp.constraints.x0 = zeros(nx,1);

% state constraints on the first, second and sixth state
% (not applied to the initial state)
ocp.constraints.idxbx = [0 1 5];                % zero-based indices
infty = get_acados_infty();                     % for one-sided constraints
ocp.constraints.lbx = [-pi/6; -pi/6; -1];       % state lower bounds
ocp.constraints.ubx = [ pi/6;  pi/6; infty];    % state upper bounds

% use the same constraints on the terminal state
ocp.constraints.idxbx_e = ocp.constraints.idxbx;
ocp.constraints.lbx_e = ocp.constraints.lbx;
ocp.constraints.ubx_e = ocp.constraints.ubx;

%% ocp options
ocp.solver_options.N_horizon = N_horizon;
ocp.solver_options.tf = tf;
ocp.solver_options.integrator_type = 'DISCRETE';

%% create the solver
ocp_solver = AcadosOcpSolver(ocp);

%% simulate the closed-loop system
x0 = zeros(nx,1);                           % initial state for the simulation
xr = [0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0];  % state reference
nsim = 15;                                  % number of simulation steps

X = nan(nx,nsim+1);             % state log
X(:,1) = x0;                    % first entry is the initial state
U = nan(nu,nsim);               % control input log
solve_time_log = nan(1,nsim);   % for solver performance evaluation

x = x0;
for isim = 1:nsim
    % set the current state
    ocp_solver.set('constr_x0', x);

    % set the reference (last argument is the stage)
    for k = 1:N_horizon-1  % intermediate stages
        ocp_solver.set('cost_y_ref', [zeros(nu,1); xr], k);
    end
    ocp_solver.set('cost_y_ref_e', xr, N_horizon);  % terminal stage

    % solve the ocp
    ocp_solver.solve();

    % check the solver output
    status = ocp_solver.get('status');
    if status ~= 0
        warning(['acados ocp solver failed with status ',num2str(status)]);
    end

    % apply the first control input to the plant, simulate one step
    ctrl = ocp_solver.get('u', 0);
    x = Ad*x + Bd*ctrl;

    % log the data
    X(:,isim+1) = x;
    U(:,isim) = ctrl;
    solve_time_log(isim) = ocp_solver.get('time_tot');
end

disp([newline,'Average solve time: ', num2str(1e3*mean(solve_time_log)), ' ms'])

%% plot the results
figure
subplot(2,1,1)
plot(X(3,:))
hold on; yline(xr(3),'r--')
ylim padded
xlim tight
ylabel('$x_3$')
subplot(2,1,2)
stairs(U')
yline([ocp.constraints.lbu(1) ocp.constraints.ubu(1)],'r--')
ylim padded
xlim tight
ylabel('$u$')
xlabel('step')

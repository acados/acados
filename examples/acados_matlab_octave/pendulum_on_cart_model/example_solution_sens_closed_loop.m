%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



clear all; clc;

model_path = fullfile(pwd,'..','pendulum_on_cart_model');
addpath(model_path)

% initial state
xcurrent = [0.2; 0; 0; 0];

%% discretization
N = 20; % number of shooting intervals
% nonuniform discretization
T = 1.0;
shooting_nodes = linspace(0, T, N+1);
h = T/N; % sampling time = length of first shooting interval

nlp_solver = 'sqp'; % sqp, sqp_rti
qp_solver = 'partial_condensing_hpipm';
% full_condensing_hpipm, partial_condensing_hpipm, full_condensing_qpoases
qp_solver_cond_N = 5; % for partial condensing

% we add some model-plant mismatch by choosing different integration
% methods for model (within the OCP) and plant:

% integrator model
model_sim_method = 'erk';
model_sim_method_num_stages = 1;
model_sim_method_num_steps = 2;

% integrator plant
plant_sim_method = 'irk';
plant_sim_method_num_stages = 3;
plant_sim_method_num_steps = 3;

%% model dynamics
model = pendulum_on_cart_model();
nx = model.nx;
nu = model.nu;

model_name = 'pendulum';

%% OCP formulation
ocp = AcadosOcp();
ocp.model.name = model_name;
ocp.model.x = model.sym_x;
ocp.model.u = model.sym_u;
ocp.model.xdot = model.sym_xdot;
ocp.model.f_expl_expr = model.dyn_expr_f_expl;
ocp.model.cost_expr_ext_cost = model.cost_expr_ext_cost;
ocp.model.cost_expr_ext_cost_e = model.cost_expr_ext_cost_e;
ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';
ocp.model.con_h_expr_0 = model.constr_expr_h;
ocp.model.con_h_expr = model.constr_expr_h;
U_max = 80;
ocp.constraints.lh_0 = -U_max;
ocp.constraints.uh_0 = U_max;
ocp.constraints.lh = -U_max;
ocp.constraints.uh = U_max;
ocp.constraints.x0 = xcurrent;

ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T;
ocp.solver_options.shooting_nodes = shooting_nodes;
ocp.solver_options.nlp_solver_type = upper(nlp_solver);
ocp.solver_options.integrator_type = upper(model_sim_method);
ocp.solver_options.sim_method_num_stages = model_sim_method_num_stages;
ocp.solver_options.sim_method_num_steps = model_sim_method_num_steps;
ocp.solver_options.qp_solver = upper(qp_solver);
ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
ocp_solver = AcadosOcpSolver(ocp);

x_traj_init = zeros(nx, N+1);
u_traj_init = zeros(nu, N);


%% plant: create integrator
sim = AcadosSim();
sim.model.name = [model_name, '_plant'];
sim.model.x = model.sym_x;
sim.model.u = model.sym_u;
sim.model.xdot = model.sym_xdot;
sim.model.f_impl_expr = model.dyn_expr_f_impl;
sim.solver_options.Tsim = h;
sim.solver_options.integrator_type = upper(plant_sim_method);
sim.solver_options.num_stages = plant_sim_method_num_stages;
sim.solver_options.num_steps = plant_sim_method_num_steps;
sim_solver = AcadosSimSolver(sim);

%% Simulation
N_sim = 100;

x_sim = zeros(nx, N_sim+1);
u_sim = zeros(nu, N_sim);

x_sim(:,1) = xcurrent;

for i=1:N_sim
    % update initial state
    xcurrent = x_sim(:,i);
    ocp_solver.set('constr_x0', xcurrent);

    if i == 1 || i == floor(N_sim/2)
        % solve
        ocp_solver.solve();
        % get solution
        u0 = ocp_solver.get('u', 0);
        status = ocp_solver.get('status'); % 0 - success
        x_lin = xcurrent;
        u_lin = u0;

        sens_u = zeros(nx, N);
        % get sensitivities w.r.t. initial state value with index
        field = 'ex'; % equality constraint on states
        stage = 0;
        for index = 0:nx-1
            ocp_solver.eval_param_sens(field, stage, index);
            sens_u(index+1,:) = ocp_solver.get('sens_u');
        end
    else
        % use feedback policy
        delta_x = xcurrent-x_lin;
        u0 = u_lin + sens_u(:, 1)' * delta_x;
    end

    % set initial state
    sim_solver.set('x', xcurrent);
    sim_solver.set('u', u0);

    % solve
    sim_solver.solve();

    % get simulated state
    x_sim(:,i+1) = sim_solver.get('xn');
    u_sim(:,i) = u0;
end

disp('final state')
format long e
disp(x_sim(:,N_sim+1))
%      1.392073955008204e-03
%      4.247720422461933e-05
%     -7.518679517751918e-05
%     -1.671811407900214e-04

%% Plots
ts = linspace(0, N_sim*h, N_sim+1);
figure; hold on;
States = {'p', 'theta', 'v', 'dtheta'};
v_mean = 0;
p_ref = ts*v_mean;

y_ref = zeros(nx, N_sim+1);
y_ref(1, :) = p_ref;
y_ref(3, :) = v_mean;

for i=1:length(States)
    subplot(length(States), 1, i);
    grid on; hold on;
    plot(ts, x_sim(i,:));
    plot(ts, y_ref(i, :));
    ylabel(States{i});
    xlabel('t [s]')
    legend('closed-loop', 'reference')
end

figure
stairs(ts(1:end), [u_sim'; u_sim(end)])
ylabel('F [N]')
xlabel('t [s]')
grid on

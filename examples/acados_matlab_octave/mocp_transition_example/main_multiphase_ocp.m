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


settings = get_example_settings();

N_list = [10, 1, 15];
n_phases = length(N_list);
N_horizon = sum(N_list);

% create_multiphase_ocp_solver
ocp = AcadosMultiphaseOcp(N_list);

phase_1 = formulate_double_integrator_ocp(settings);
ocp.set_phase(phase_1, 1);

phase_2 = AcadosOcp();
phase_2.model = get_transition_model();

% define transition cost
phase_2.cost.cost_type = 'NONLINEAR_LS';
phase_2.model.cost_y_expr = phase_2.model.x;
phase_2.cost.W = diag([settings.L2_COST_P, 1e-1 * settings.L2_COST_V]);
phase_2.cost.yref = zeros(2, 1);
ocp.set_phase(phase_2, 2);

phase_3 = formulate_single_integrator_ocp(settings);
ocp.set_phase(phase_3, 3);

% set mocp specific options
ocp.mocp_opts.integrator_type = {'ERK', 'DISCRETE', 'ERK'};

% set solver options, common for AcadosOcp and AcadosMultiphaseOcp
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.tf = settings.T_HORIZON;
T_HORIZON_1 = 0.4 * settings.T_HORIZON;
T_HORIZON_2 = settings.T_HORIZON - T_HORIZON_1;
ocp.solver_options.time_steps = [T_HORIZON_1 / N_list(1) * ones(1, N_list(1)), ...
                                0.0, ...  % transition stage
                                T_HORIZON_2 / N_list(3) * ones(1, N_list(3))];

ocp.solver_options.cost_scaling = [T_HORIZON_1 / N_list(1) * ones(1, N_list(1)), ...
                                1.0, ...  % transition stage
                                T_HORIZON_2 / N_list(3) * ones(1, N_list(3)), ...
                                1.0];  % terminal cost

ocp.solver_options.store_iterates = true;

ocp_solver = AcadosOcpSolver(ocp);

% initialize x trajectory using flattened format
x0 = ocp.constraints{1}.x0;
x_init = [repmat(x0, 1, N_list(1)+N_list(2)) repmat(x0(1), 1, N_list(3)+1)];
ocp_solver.set('x', x_init);

% update state bounds using flattened format
lbx = [repmat([-10, -5], 1, N_list(1)) repmat([-10], 1, N_list(3)+1)];
ocp_solver.set('constr_lbx', lbx);

% need to set initial state after updating the bounds on x as this again overwrites lbx_0
ocp_solver.set('constr_x0', x0);

ocp_solver.solve();
ocp_solver.print();

iterate = ocp_solver.get_iterate(ocp_solver.get('sqp_iter'));

%% extract solution
x_traj = cell(N_horizon+1, 1);
u_traj = cell(N_horizon, 1);
for i=0:N_horizon
    x_traj{i+1, 1} = ocp_solver.get('x', i);
end
for i=0:N_horizon-1
    u_traj{i+1, 1} = ocp_solver.get('u', i);
end

if iterate.x_traj{end} ~= x_traj{end}
    error("x and last iterate for x should be the same.")
end

t_grid = ocp.solver_options.shooting_nodes;

p_traj_1 = zeros(N_list(1)+1, 1);
v_traj_1 = zeros(N_list(1)+1, 1);
for i=1:N_list(1)+1
    p_traj_1(i) = x_traj{i}(1);
    v_traj_1(i) = x_traj{i}(2);
end

p_traj_3 = zeros(N_list(3)+1, 1);
for i=1:N_list(3)+1
    p_traj_3(i) = x_traj{N_list(1)+N_list(2)+i}(1);
end

v_traj_3 = zeros(N_list(3), 1);
for i=1:N_list(3)
    v_traj_3(i) = u_traj{N_list(1)+N_list(2)+i}(1);
end


t_grid_phases = cell(n_phases, 1);
for i=1:n_phases
    t_grid_phases{i} = t_grid(sum(ocp.start_idx(i))+1:ocp.end_idx(i)+1);
end

a_traj = zeros(N_list(1), 1);
for i=1:N_list(1)
    a_traj(i) = u_traj{i}(1);
end

%% plot trajectories in subplots
figure;

subplot(3, 1, 1);
hold on;
plot(t_grid_phases{1}, p_traj_1(:, 1), '-');
plot(t_grid_phases{3}, p_traj_3(:, 1), '-');
legend('phase 0', 'phase 2')
ylabel('position');
xlim([0, settings.T_HORIZON]);

subplot(3, 1, 2);
hold on;
plot(t_grid_phases{1}, v_traj_1(:, 1), '-');
stairs(t_grid_phases{3}, [v_traj_3; v_traj_3(end)], '-');
legend('phase 0', 'phase 2')
ylabel('velocity');
xlim([0, settings.T_HORIZON]);

subplot(3, 1, 3);
stairs(t_grid_phases{1}, [a_traj; a_traj(end)], '-');
ylabel('acceleration');
xlabel('t [s]');
xlim([0, settings.T_HORIZON]);

% to show plot in octave, commented to run example on CI
% if is_octave()
%     waitforbuttonpress;
% end

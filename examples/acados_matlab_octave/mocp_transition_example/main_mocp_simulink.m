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
settings.WITH_X_BOUNDS = false;

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
% add dummy constraints to test Simulink
phase_3.model.con_h_expr = [phase_3.model.x; phase_3.model.x^2];
phase_3.constraints.lh = [-100; -200];
phase_3.constraints.uh = [100; 200];

phase_3.constraints.idxbx = [0];
phase_3.constraints.lbx = -200;
phase_3.constraints.ubx = 200;

% values for Simulink
N_3 = 15;
lbx_all = repmat(phase_3.constraints.lbx, N_3);
ubx_all = repmat(phase_3.constraints.ubx, N_3);
lh_all = repmat(phase_3.constraints.lh, N_3);
uh_all = repmat(phase_3.constraints.uh, N_3);

ocp.set_phase(phase_3, 3);

% set mocp specific options
ocp.mocp_opts.integrator_type = {'ERK', 'DISCRETE', 'ERK'};

% set solver options, common for AcadosOcp and AcadosMultiphaseOcp
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.tf = settings.T_HORIZON + 1.0;
T_HORIZON_1 = 0.4 * settings.T_HORIZON;
T_HORIZON_2 = settings.T_HORIZON - T_HORIZON_1;
ocp.solver_options.time_steps = [T_HORIZON_1 / N_list(1) * ones(1, N_list(1)), ...
                                1.0, ...  % transition stage
                                T_HORIZON_2 / N_list(3) * ones(1, N_list(3))];

%% simulink options
simulink_opts = get_acados_simulink_opts_mocp();
simulink_opts.inputs.x_init = 1;
simulink_opts.inputs.u_init = 1;
simulink_opts.inputs.pi_init = 1;
simulink_opts.inputs.ignore_inits = 1;
simulink_opts.inputs.reset_solver = 1;
simulink_opts.inputs.lbx = 1;
simulink_opts.inputs.ubx = 1;
simulink_opts.inputs.lh = 1;
simulink_opts.inputs.uh = 1;

simulink_opts.outputs.xtraj = 1;
simulink_opts.outputs.utraj = 1;
simulink_opts.outputs.pi_all = 1;
simulink_opts.outputs.CPU_time = 0;

ocp.simulink_opts = simulink_opts;

x0_original = phase_1.constraints.x0;

tol = 1e-12;
ocp.solver_options.nlp_solver_tol_stat = tol;
ocp.solver_options.nlp_solver_tol_eq = tol;
ocp.solver_options.nlp_solver_tol_ineq = tol;
ocp.solver_options.nlp_solver_tol_comp = tol;

%% create solver
ocp_solver = AcadosOcpSolver(ocp);

%% reference solve
ocp_solver.solve();
ocp_solver.print()


%% extract solution
x_traj = [];
u_traj = [];
pi_traj = [];
for i=0:N_horizon
    x_traj = [x_traj; ocp_solver.get('x', i)];
end
for i=0:N_horizon-1
    u_traj = [u_traj; ocp_solver.get('u', i)];
    pi_traj = [pi_traj; ocp_solver.get('pi', i)];
end
%% simulink test
cd c_generated_code
make_sfun; % ocp solver
cd ..;
n_sim = 3;

%% Test Simulink example block
for itest = 1:2
    % open_system('mocp_simulink_block')
    lbu_all = [repmat(phase_1.constraints.lbu, N_list(1), 1);
               repmat(phase_2.constraints.lbu, N_list(2), 1);
               repmat(phase_3.constraints.lbu, N_list(3), 1)];
    ubu_all = [repmat(phase_1.constraints.ubu, N_list(1), 1);
               repmat(phase_2.constraints.ubu, N_list(2), 1);
               repmat(phase_3.constraints.ubu, N_list(3), 1)];

    lbx_all = [repmat(phase_1.constraints.lbx, N_list(1)-1, 1);
               repmat(phase_2.constraints.lbx, N_list(2), 1);
               repmat(phase_3.constraints.lbx, N_list(3), 1)];
    ubx_all = [repmat(phase_1.constraints.ubx, N_list(1)-1, 1);
               repmat(phase_2.constraints.ubx, N_list(2), 1);
               repmat(phase_3.constraints.ubx, N_list(3), 1)];

    lh_all = [repmat(phase_1.constraints.lh, N_list(1)-1, 1);
               repmat(phase_2.constraints.lh, N_list(2), 1);
               repmat(phase_3.constraints.lh, N_list(3), 1)];
    uh_all = [repmat(phase_1.constraints.uh, N_list(1)-1, 1);
               repmat(phase_2.constraints.uh, N_list(2), 1);
               repmat(phase_3.constraints.uh, N_list(3), 1)];


    bu_fact = 1;
    if itest == 2
        bu_fact = 0.8;
        lbu_all = bu_fact * lbu_all;
        ubu_all = bu_fact * ubu_all;
    end

    for stage = 0:N_list(1)-1
        ocp_solver.set('constr_lbu', bu_fact * phase_1.constraints.lbu, stage);
        ocp_solver.set('constr_ubu', bu_fact * phase_1.constraints.ubu, stage);
    end
    for stage = N_list(1):(N_list(1)+N_list(2)-1)
        ocp_solver.set('constr_lbu', bu_fact * phase_2.constraints.lbu, stage);
        ocp_solver.set('constr_ubu', bu_fact * phase_2.constraints.ubu, stage);
    end
    for stage = (N_list(1)+N_list(2)):(sum(N_list)-1)
        ocp_solver.set('constr_lbu', bu_fact * phase_3.constraints.lbu, stage);
        ocp_solver.set('constr_ubu', bu_fact * phase_3.constraints.ubu, stage);
    end
    % reference solve
    ocp_solver.solve();
    % extract solution
    x_traj = [];
    u_traj = [];
    pi_traj = [];
    for i=0:N_horizon
        x_traj = [x_traj; ocp_solver.get('x', i)];
    end
    for i=0:N_horizon-1
        u_traj = [u_traj; ocp_solver.get('u', i)];
        pi_traj = [pi_traj; ocp_solver.get('pi', i)];
    end

    x_init = 0 * x_traj;
    u_init = 0 * u_traj;
    pi_init = 0 * pi_traj;
    ignore_inits = 0;
    reset_solver = 0;

    %
    out_sim = sim('mocp_simulink_block', 'SaveOutput', 'on');
    fprintf('\nSuccessfully ran simulink block.\n');

    % Evaluation
    fprintf('\nTest results on SIMULINK simulation.\n')

    %
    disp('checking KKT residual, should be zero just because not evaluated')
    kkt_signal = out_sim.logsout.getElement('KKT_residual');
    kkt_val = kkt_signal.Values.data;
    if any(kkt_val > 1e-6)
        error('failed');
    end

    sqp_iter_signal = out_sim.logsout.getElement('sqp_iter');
    sqp_iter_val = sqp_iter_signal.Values.Data;
    disp('checking sqp_iter_val, should be 1 because problem is a QP.')
    if any(sqp_iter_val ~= 1)
        error('failed');
    end

    status_signal = out_sim.logsout.getElement('status');
    status_val = status_signal.Values.Data;
    disp('checking status, should be 0 (success).')
    if any(status_val ~= 0)
        error(['failed. got status values:' mat2str(status_val)]);
    end

    utraj_signal = out_sim.logsout.getElement('utraj');
    u_simulink = utraj_signal.Values.Data(1, :);
    disp('checking u values.')
    if norm(u_simulink(:) - u_traj(:)) > 1e-8
        error('failed');
    end

    xtraj_signal = out_sim.logsout.getElement('xtraj');
    x_simulink = xtraj_signal.Values.Data(1, :);
    disp('checking x values, should match solution in MATLAB.')
    if norm(x_simulink(:) - x_traj(:)) > 1e-8
        error('failed');
    end
    %
    pi_signal = out_sim.logsout.getElement('pi_all');
    pi_simulink = pi_signal.Values.Data(1, :);
    disp('checking pi values, should match solution in MATLAB.')
    if norm(pi_simulink(:) - pi_traj(:)) > 1e-8
        error('failed');
    end

    x1_signal = out_sim.logsout.getElement('x1');
    x1_val = x1_signal.Values.data(1, :);
    x1_ref = ocp_solver.get('x', 1);
    disp('checking x1_value')
    if norm(x1_val(:) - x1_ref(:)) > 1e-8
        error('failed');
    end

    u0_signal = out_sim.logsout.getElement('u0');
    u0_val = u0_signal.Values.data(1, :);
    u0_ref = ocp_solver.get('u', 0);
    disp('checking u0_value')
    if norm(u0_val(:) - u0_ref(:)) > 1e-8
        error('failed');
    end

end

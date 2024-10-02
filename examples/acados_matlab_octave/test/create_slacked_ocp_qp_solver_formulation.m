function [ocp, simulink_opts] = create_slacked_ocp_qp_solver_formulation(N)
import casadi.*

%% solver settings
T = 1; % [s] prediction horizon length
nlp_solver = 'sqp'; % sqp, sqp_rti
% integrator type
sim_method = 'erk'; % erk, irk, irk_gnsf

nx = 3; % state size
nu = 3; % input size
ny = nx+nu;
ny_e = nx;

x0 = ones(nx, 1); % initial state

% OCP formulation object
ocp = AcadosOcp();
%% model
ocp.model.x = SX.sym('x', nx);
ocp.model.u = SX.sym('u', nu);
ocp.model.name = 'soft_constraints_qp';
ocp.model.f_expl_expr = ocp.model.u;

% cost
W = eye(ny);
Vx = zeros(ny,nx); Vx(1:nx,:) = eye(nx);        % state-to-output matrix in lagrange term
Vu = zeros(ny,nu); Vu(nx+1:ny,:) = eye(nu);     % input-to-output matrix in lagrange term
Vx_e = zeros(ny_e,nx); Vx_e(1:nx,:) = eye(nx);  % state-to-output matrix in mayer term

ocp.cost.cost_type = 'LINEAR_LS';
ocp.cost.cost_type_e = 'LINEAR_LS';
ocp.cost.W = W;
ocp.cost.Vx = Vx;
ocp.cost.Vu = Vu;
ocp.cost.Vx_e = Vx_e;

ocp.cost.W_e = 5 * W(1:ny_e,1:ny_e);
ocp.cost.yref = zeros(ny,1);
ocp.cost.yref_e = zeros(ny_e,1);


% constraints
ocp.constraints.x0 = x0;
% soft constraints
ocp.constraints.lbu = -0.5*ones(nu,1);
ocp.constraints.ubu = 0.5*ones(nu,1);
ocp.constraints.idxbu = (1:nu) - 1;
ocp.constraints.idxsbu = (1:nu) - 1;

ns = nu;
slack_penalty = 1e2;
ocp.cost.Zl_0 = slack_penalty*ones(ns,1);
ocp.cost.Zu_0 = slack_penalty*ones(ns,1);
ocp.cost.zl_0 = slack_penalty*ones(ns,1);
ocp.cost.zu_0 = slack_penalty*ones(ns,1);

ocp.cost.Zl = slack_penalty*ones(ns,1);
ocp.cost.Zu = slack_penalty*ones(ns,1);
ocp.cost.zl = slack_penalty*ones(ns,1);
ocp.cost.zu = slack_penalty*ones(ns,1);

% options
ocp.solver_options.N_horizon = N;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.tf = T;
ocp.solver_options.integrator_type = 'ERK';

%% Simulink opts
simulink_opts = get_acados_simulink_opts;

% inputs
simulink_opts.inputs.y_ref = 0;
simulink_opts.inputs.y_ref_0 = 0;
simulink_opts.inputs.y_ref_e = 0;

simulink_opts.inputs.cost_W = 0;
simulink_opts.inputs.cost_W_0 = 0;
simulink_opts.inputs.cost_W_e = 0;

simulink_opts.inputs.lbx = 0;
simulink_opts.inputs.ubx = 0;
simulink_opts.inputs.lbx_e = 0;
simulink_opts.inputs.ubx_e = 0;
simulink_opts.inputs.lbu = 0;
simulink_opts.inputs.ubu = 0;

% soft constraints ports
simulink_opts.inputs.cost_zl = 1;
simulink_opts.inputs.cost_zu = 1;
simulink_opts.inputs.cost_Zl = 1;
simulink_opts.inputs.cost_Zu = 1;

simulink_opts.inputs.slacks_init = 1;


simulink_opts.outputs.slack_values = 1;

% outputs
simulink_opts.outputs.u0 = 0;
simulink_opts.outputs.utraj = 1;
simulink_opts.outputs.xtraj = 1;
simulink_opts.outputs.pi_all = 0;

simulink_opts.outputs.solver_status = 1;
simulink_opts.outputs.CPU_time = 0;
simulink_opts.outputs.sqp_iter = 1;
simulink_opts.outputs.x1 = 0;
end
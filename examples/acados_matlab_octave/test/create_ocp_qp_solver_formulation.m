function [ocp_model, ocp_opts, simulink_opts, x0] = create_ocp_qp_solver_formulation(N)
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

%% acados ocp model
ocp_model = acados_ocp_model();

sym_x = SX.sym('x', nx);
sym_u = SX.sym('u', nu);


ocp_model.set('name', 'trivial');
ocp_model.set('T', T);  % prediction horizon

% symbolics
ocp_model.set('sym_x', sym_x);
ocp_model.set('sym_u', sym_u);

% cost
cost_type = 'linear_ls';
cost_type_e = 'linear_ls';
Vx = zeros(ny,nx); Vx(1:nx,:) = eye(nx);        % state-to-output matrix in lagrange term
Vu = zeros(ny,nu); Vu(nx+1:ny,:) = eye(nu);     % input-to-output matrix in lagrange term
Vx_e = zeros(ny_e,nx); Vx_e(1:nx,:) = eye(nx);  % state-to-output matrix in mayer term
W = eye(ny);
W_e = 5 * W(1:ny_e,1:ny_e);                         % cost weights in mayer term
y_ref = zeros(ny,1);                            % set intial references
y_ref_e = zeros(ny_e,1);

% cost
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type_e);
ocp_model.set('cost_Vx', Vx);
ocp_model.set('cost_Vu', Vu);
ocp_model.set('cost_Vx_e', Vx_e);
ocp_model.set('cost_W', W);
ocp_model.set('cost_W_e', W_e);
ocp_model.set('cost_y_ref', y_ref);
ocp_model.set('cost_y_ref_e', y_ref_e);

% dynamics
if (strcmp(sim_method, 'erk'))
    ocp_model.set('dyn_type', 'explicit');
    ocp_model.set('dyn_expr_f', sym_u);
end

% constraints
ocp_model.set('constr_x0', x0);  % set the initial state

%% acados ocp options
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', nlp_solver);
ocp_opts.set('sim_method', sim_method);


%% Simulink opts
simulink_opts = get_acados_simulink_opts;

% inputs
simulink_opts.inputs.y_ref = 0;
simulink_opts.inputs.y_ref_0 = 0;
simulink_opts.inputs.y_ref_e = 0;

simulink_opts.inputs.cost_W = 0;
simulink_opts.inputs.cost_W_0 = 0;
simulink_opts.inputs.cost_W_e = 0;

simulink_opts.inputs.x_init = 1;
simulink_opts.inputs.u_init = 1;
simulink_opts.inputs.pi_init = 1;

simulink_opts.inputs.reset_solver = 1;
simulink_opts.inputs.ignore_inits = 1;

% outputs
simulink_opts.outputs.u0 = 0;
simulink_opts.outputs.utraj = 1;
simulink_opts.outputs.xtraj = 1;
simulink_opts.outputs.pi_all = 1;

simulink_opts.outputs.solver_status = 1;
simulink_opts.outputs.CPU_time = 0;
simulink_opts.outputs.sqp_iter = 1;
simulink_opts.outputs.x1 = 0;
end
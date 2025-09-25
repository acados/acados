import casadi.*

%% Circuit Paramters - Per Unit Values
Rpu_load = 1;
Cpu_out = 97.9;
Lpu_out = 0.562;
Vpu_source = 1;
Tpu_sampling = 1;

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% system model
model_name = 'buckConv_dae';

% Differential States (x vector)
x1 = SX.sym('x1');  % inductor current
x2 = SX.sym('x2');  % capacitor voltage output
x = vertcat(x1, x2); % state vector

% algebraic states
z = SX.sym('z');

% controls
u = SX.sym('u');  % duty cycle input

% state derivatives (x dot vector)
x1_dot = SX.sym('x1_dot'); % derivative of inductor current
x2_dot = SX.sym('x2_dot'); % derivative of capacitor voltage
x_dot = vertcat(x1_dot, x2_dot); % state derivative vector

%% Dynamics: implicit DAE formulation (index-1)
f_impl = vertcat(- x1_dot + u * (Vpu_source/Lpu_out) - x2 * (1/Lpu_out), ...
                 x1 * (1/Cpu_out) - x2 * (1/(Rpu_load*Cpu_out)) - x2_dot,...
                 z - (u * x1));

nx = length(x);
nu = length(u);
nz = length(z);

ny = nu+nx+nz; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term

% initialize values cost selection matrices
Vx = zeros(ny,nx);
Vz = zeros(ny,nz);
Vu = zeros(ny,nu);

Vx(2,2) = 1; % select none state to track in cost function

% terminal cost (mayer) terms (no inputs or algebraic states)
Vx_e = Vx(1:nx,:);

% now make the reference vector for that reference tracking
iL_ref = 0.1;
y_ref = zeros(ny,1); % initialize all values to 0
y_ref(2) = iL_ref; % set reference to desired value

y_ref_e = y_ref(1:nx); % terminal reference, now same as path but without input/algebraic reference (because no input based terms in terminal cost)

% Weight matrix for reference tracking
weights = zeros(ny,1); % initialize correct length weight vector
weights(2) = 1; % set weights we want to track
W = diag(weights); % make square diagonal matrix from vector

weights_e = weights(1:nx); % terminal weights only on states
W_e = diag(weights_e);


model = AcadosModel();
model.x = x;
model.xdot = x_dot;
model.u = u;
model.z = z;

model.f_impl_expr = f_impl;
model.name = model_name;

%% OCP formulation object
ocp = AcadosOcp();
ocp.model = model;

ocp.cost.cost_type = 'EXTERNAL'; % linear_ls, ext_cost
ocp.cost.cost_type_e = 'EXTERNAL'; % linear_ls, ext_cost

y = Vx * x + Vu * u + Vz * z; % output vector
y_e = Vx_e * x; % terminal output vector

ocp.model.cost_expr_ext_cost = 0.5* (y - y_ref)'*W*(y - y_ref); % external cost expression
ocp.model.cost_expr_ext_cost_e = 0.5* (y_e - y_ref_e)'*W_e*(y_e - y_ref_e); % external cost expression

% constraints on input
nbu = nu; % number of bounds on controls u
idxbu = 0:(nu-1); % which controls the bounds are on
lbu = zeros(nu, 1); % lower bound for constraint
ubu = ones(nu, 1);  % upper bound for constraint
ocp.constraints.lbu = lbu;
ocp.constraints.ubu = ubu;
ocp.constraints.idxbu = idxbu;

% Combined Constraint on state/input
const_lh = -1; % current constraint in pu, lower bound
const_uh = 5; % current constraint in pu, upper bound

x0 = zeros(nx,1); % initial condition
ocp.constraints.x0 = x0;
ocp.model.con_h_expr = x1;
ocp.constraints.lh = const_lh;
ocp.constraints.uh = const_uh;

ocp_N = 80;
T = ocp_N*Tpu_sampling; % Prediction horizon in time form for acados
ocp.solver_options.tf = T;
ocp.solver_options.N_horizon = ocp_N;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
ocp.solver_options.integrator_type = 'IRK';
ocp.solver_options.sim_method_num_stages = 2;
ocp.solver_options.sim_method_num_steps = 1;
ocp.solver_options.sim_method_newton_tol = 1e-6;

% create solver
ocp_solver = AcadosOcpSolver(ocp);

x_SS = [iL_ref;iL_ref*Rpu_load];
u_SS = iL_ref*Rpu_load/Vpu_source;
z_SS = iL_ref*u_SS;

x_traj_init = repmat(x0, 1, ocp_N + 1);
u_traj_init = u_SS*ones(nu, ocp_N);
z_traj_init = z_SS*ones(nz, ocp_N);
xdot_traj_init = 0*ones(nx, ocp_N);

% set initial trajectories
ocp_solver.set('init_x', x_traj_init);
ocp_solver.set('init_u', u_traj_init);
ocp_solver.set('init_z', z_traj_init); % only set at initial time
ocp_solver.set('init_xdot', xdot_traj_init); % only set at initial time

ocp_solver.solve();
ocp_solver.print('stat')

% check status
status = ocp_solver.get('status');
sqp_iter = ocp_solver.get('sqp_iter');
qp_iter = ocp_solver.get('qp_iter');
time_tot = ocp_solver.get('time_tot');
time_lin = ocp_solver.get('time_lin');
time_qp_sol = ocp_solver.get('time_qp_sol');

fprintf('\nstatus = %d, sqp_iter = %d, qp_iter = %d, time_int = %f [ms] (time_lin = %f [ms], time_qp_sol = %f [ms])\n',...
    status, sqp_iter, qp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3);

if status~=0
    disp('acados ocp_solver solver failed');
    keyboard
end

x_traj = ocp_solver.get('x');
u_traj = ocp_solver.get('u');
z_traj = ocp_solver.get('z');

xu_product = x_traj(1, 1:end-1) .* u_traj;

max_violation = max(xu_product - z_traj)
disp(['Max violation of algebraic constraint' max_violation])

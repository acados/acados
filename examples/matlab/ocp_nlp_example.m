import casadi.*
import acados.*

N = 10;
[ode_fun, nx, nu] = chen_model();
nlp = ocp_nlp(struct('N', N, 'nx', nx, 'nu', nu));

% ODE Model
step = 0.1;
nlp.set_model(ode_fun, step);

% Cost function
Q = diag([1.0, 1.0]);
R = 1e-2;
x_ref = [0.01; 0.2];
x = SX.sym('x',nx);
u = SX.sym('u',nu);
uN = SX.sym('u',0);
F = ocp_nlp_function(Function('ls_cost', {x,u}, {vertcat(x,u)}));
FN = ocp_nlp_function(Function('ls_costN', {x,uN}, {x}));
cost_functions = cell(N+1, 1);
cost_matrices = cell(N+1, 1);
cost_ref = cell(N+1, 1);
for i=1:N
    cost_functions{i} = F;
    cost_matrices{i} = blkdiag(Q, R);
    cost_ref{i} = [x_ref; 0];
end
cost_functions{N+1} = FN;
cost_matrices{N+1} = Q;
cost_ref{N+1} = x_ref;
ls_cost = ocp_nlp_ls_cost(N, cost_functions);
ls_cost.ls_cost_matrix = cost_matrices;
ls_cost.ls_cost_ref = cost_ref;
nlp.set_cost(ls_cost);

% Constraints
G = ocp_nlp_function(Function('path_constraints', {x, u}, {SX([])}));
GN = ocp_nlp_function(Function('path_constraintsN', {x, uN}, {SX([])}));
path_constraints = cell(N+1, 1);
for i=1:N
    path_constraints{i} = G;
end
path_constraints{N+1} = GN;
nlp.set_path_constraints(path_constraints);
 
solver = ocp_nlp_solver('sqp', nlp, struct('integrator_steps', 2, 'qp_solver', 'condensing_qpoases', 'sensitivity_method', 'gauss-newton'));
 
% Simulation
num_iters = 1000;
STATES = zeros(nx, num_iters+1);
STATES(:, 1) = [0.1; 0.1];
for i=1:num_iters
    output = solver.evaluate(STATES(:, i));
    STATES(:, i+1) = output.states{2};
end

plot(STATES(1, :), STATES(2, :));
axis equal
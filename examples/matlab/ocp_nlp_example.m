import casadi.*
import acados.*

N = 10;
[ode_fun, nx, nu] = chen_model();
ng = cell(N+1, 1);
for i=1:N
    ng{i} = 1;
end
ng{N+1} = 0;
nlp = ocp_nlp(struct('N', N, 'nx', nx, 'nu', nu, 'ng', {ng}));

% ODE Model
step = 0.1;
nlp.set_model(ode_fun, step);

% Cost function
Q = diag([1.0, 1.0]);
R = 1e-2;
x = SX.sym('x',nx);
u = SX.sym('u',nu);
u_N = SX.sym('u',0);
f = ocp_nlp_function(Function('ls_cost', {x, u}, {vertcat(x, u)}));
f_N = ocp_nlp_function(Function('ls_cost_N', {x, u_N}, {x}));
cost_functions = cell(N+1, 1);
cost_matrices = cell(N+1, 1);
for i=1:N
    cost_functions{i} = f;
    cost_matrices{i} = blkdiag(Q, R);
end
cost_functions{N+1} = f_N;
cost_matrices{N+1} = Q;

ls_cost = ocp_nlp_ls_cost(N, cost_functions);
ls_cost.ls_cost_matrix = cost_matrices;
nlp.set_cost(ls_cost);

% Constraints
g = ocp_nlp_function(Function('path_constraint', {x, u}, {u}));
g_N = ocp_nlp_function(Function('path_constraint_N', {x, u_N}, {SX([])}));
path_constraints = cell(N+1, 1);
for i=1:N
    path_constraints{i} = g;
    nlp.lg{i} = -0.5;
    nlp.ug{i} = +0.5;
end
path_constraints{N+1} = g_N;
nlp.set_path_constraints(path_constraints);
 
solver = ocp_nlp_solver('sqp', nlp, struct('integrator_steps', 2, 'qp_solver', 'condensing_qpoases', 'sensitivity_method', 'gauss-newton'));
 
% Simulation
num_iters = 20;
STATES = zeros(nx, num_iters+1);
STATES(:, 1) = [0.1; 0.1];
for i=1:num_iters
    output = solver.evaluate(STATES(:, i));
    STATES(:, i+1) = output.states{2};
end

plot(STATES(1, :), STATES(2, :));
axis equal
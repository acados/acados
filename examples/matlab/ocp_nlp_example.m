import acados.*

N = 10;
[nx, nu, ode_fun] = example_model();
nlp = ocp_nlp(struct('N', N, 'nx', nx, 'nu', nu));

% ODE Model
step = 0.1;
nlp.set_model(ode_fun, step);

% Cost function
Q = diag([1.0, 1.0]);
R = 1e-2;
cost_matrices = cell(N+1, 1);
for i=1:N
    cost_matrices{i} = blkdiag(Q, R);
end
cost_matrices{N+1} = Q;
nlp.ls_cost_matrix = cost_matrices;

solver = ocp_nlp_solver('gauss-newton-sqp', nlp, struct('integrator_steps', 2));

% Simulation
num_iters = 50;
STATES = zeros(nx, num_iters+1);
STATES(:, 1) = [0.1; 0.1];
for i=1:num_iters
    output = solver.evaluate(STATES(:, i));
    STATES(:, i+1) = output.states{2};
end

plot(STATES(1, :), STATES(2, :));
axis equal

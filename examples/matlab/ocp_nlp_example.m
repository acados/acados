import acados.*
import casadi.*

N = 10;
nx = 2;
nu = 1;

nlp = ocp_nlp(struct('N', N, 'nx', nx, 'nu', nu))
% Specify initial condition
current_state = [0.1, 0.1];
nlp.lb{1} = current_state;
nlp.ub{1} = current_state;
% Weighting matrix
Q = diag([1.0, 1.0]);
R = 1e-2;
cost_matrices = {};
for i=1:N
    cost_matrices{i} = blkdiag(Q, R);
end
cost_matrices{N+1} = Q;
nlp.ls_cost_matrix = cost_matrices;

% The following ODE model comes from Chen1998
x = SX.sym('x', nx);
u = SX.sym('u', nu);

mu = 0.5;
rhs = vertcat(x[1] + u*(mu + (1.-mu)*x[0]), x[0] + u*(mu - 4.*(1.-mu)*x[1]));

Sx = SX.sym('Sx', nx, nx);
Su = SX.sym('Su', nx, nu);

vde_x = jtimes(rhs, x, Sx);
vde_u = jacobian(rhs, u) + jtimes(rhs, x, Su);
vde = Function('vde', [x, Sx, Su, u], [rhs, vde_x, vde_u]);
vde.generate('vde.c');
nlp.set_model('vde');

solver = ocp_nlp_solver('gauss-newton-sqp', nlp);

STATES = current_state.';

for i=1:11:
    optimal_states = solver.solve(current_state);
    current_state = optimal_states{2};
    STATES = [STATES; current_state.'];

plot(STATES(:,1), STATES(:,2));


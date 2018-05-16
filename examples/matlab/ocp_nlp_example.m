
import casadi.*
import acados.*

[ode_fun, nx, nu] = chen_model();

N = 15;

nlp = ocp_nlp(N, nx, nu);
nlp.set_dynamics(ode_fun, struct('integrator', 'rk4', 'step', 0.1));

q = 0.5; r = 1;
P = [16.5926, 11.5926; 11.5926, 16.5926];

x = SX.sym('x', nx);
u = SX.sym('u', nu);
res = Function('r', {x, u}, {vertcat(x, u)});

nlp.set_stage_cost(eye(nx+nu), zeros(nx+nu), diag([q, q, r]));
nlp.set_terminal_cost(eye(nx), zeros(nx, 1), P);

x0 = [1; 1];
nlp.set_field('lbx', 0, x0);
nlp.set_field('ubx', 0, x0);

nlp.initialize_solver('sqp', struct('qp_solver', 'hpipm'));

output = nlp.solve();

disp('states:')
disp(output.states());

disp('controls:')
disp(output.controls());

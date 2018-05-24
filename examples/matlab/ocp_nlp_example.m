
import casadi.*
import acados.*

%% Model

nx = 2;
nu = 1;

x = SX.sym('x', nx);
u = SX.sym('u', nu);
ode_fun = Function('ode_fun', {x, u}, {vertcat(x(2), u)});

N = 15;

%% NLP solver

nlp = ocp_nlp(N, nx, nu);
nlp.set_dynamics(ode_fun, struct('integrator', 'rk4', 'step', 0.1));

q = 1; r = 1;
P = eye(nx);

nlp.set_stage_cost(eye(nx+nu), zeros(nx+nu), diag([q, q, r]));
nlp.set_terminal_cost(eye(nx), zeros(nx, 1), P);

x0 = [1; 1];
nlp.set_field('lbx', 0, x0);
nlp.set_field('ubx', 0, x0);

nlp.initialize_solver('sqp', struct('qp_solver', 'hpipm'));

output = nlp.solve();

%% Plotting

states = output.states();
X = [states{1}.'];
for i=1:numel(states)-1
    X = [X; states{i+1}.'];
end

controls = output.controls();
U = [controls{1}.'];
for i=1:numel(controls)-1
    U = [U; controls{i+1}.'];
end

figure(1); clf;
subplot(211)
plot(X);
subplot(212)
plot(U);

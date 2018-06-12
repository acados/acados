
import casadi.*
import acados.*

%% Model

[ode_fun, nx, nu] = pendulum_model();

N = 20;
Ts = 0.05;
Q = diag([1, 1, 1e-10, 1e-10]);
W = blkdiag(Q, 1e-2);
WN = 1000*Q;

%% NLP solver

nlp = ocp_nlp(N, nx, nu);
nlp.set_dynamics(ode_fun, struct('integrator', 'rk4', 'step', Ts));

nlp.set_stage_cost(eye(nx+nu), zeros(nx+nu), W);
nlp.set_terminal_cost(eye(nx), zeros(nx, 1), WN);

x0 = [0; pi; 0; 0];
nlp.set_field('lbx', 0, x0);
nlp.set_field('ubx', 0, x0);

nlp.set_field('lbu', -8);
nlp.set_field('ubu', +8);

nlp.initialize_solver('sqp', struct('qp_solver', 'qpoases'));

output = nlp.solve(x0, 0);

for i=1:99

    states = output.states();
    disp(states{1});

    nlp.set_field('lbx', 0, states{2});
    nlp.set_field('ubx', 0, states{2});

    output = nlp.solve();
end

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

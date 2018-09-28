
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

nlp.initialize_solver('rti', struct('qp_solver', 'hpipm'));

%% Simulation
x = SX.sym('x', nx); u = SX.sym('u', nu);
pendulum = integrator('pendulum', 'cvodes', struct('x', x, 'p', u, 'ode', ode_fun(x, u)), struct('tf', Ts));
sim_states = x0.'; sim_controls = [];

output = nlp.solve(x0, 0);

for i=1:99

    controls = output.controls();
    sim_controls = [sim_controls; controls{1}.'];
    
    integrator_out = pendulum('x0', sim_states(end, :).', 'p', controls{1});
    sim_states = [sim_states; full(integrator_out.xf).'];
    
    nlp.set_field('lbx', 0, sim_states(end, :).');
    nlp.set_field('ubx', 0, sim_states(end, :).');

    output = nlp.solve();
end

%% Plotting
 
figure(1); clf;
subplot(211)
plot(sim_states);
subplot(212)
plot(sim_controls);

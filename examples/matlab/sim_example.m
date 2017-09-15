import acados.*

[ode_fun, nx, ~] = pendulum_model();

dt = 0.2;
sim_options = struct('time_step', dt, 'order', 4);

integrator = sim_solver('erk', ode_fun, sim_options);

current_state = [0.0; pi; 0.0; 0.0]; % pendulum hangs down
simulation = current_state;

time_grid = 0:dt:10;
for t=time_grid(1:end-1)
    % apply periodic force to pendulum with period 5s.
    output = integrator.evaluate(current_state, cos(2*pi*t/5));
    current_state = output.final_state;
    simulation = [simulation, current_state];
end

plot(time_grid, simulation.');
legend('p', 'theta', 'v', 'omega');
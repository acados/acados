settings = get_example_settings();

model = get_double_integrator_model();

ocp = formulate_double_integrator_ocp(settings);


ocp.solver_options.tf = settings.T_HORIZON;
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.N_horizon = settings.N_HORIZON;

ocp_solver = AcadosOcpSolver(ocp);

ocp_solver.solve();
ocp_solver.print()


x_traj = zeros(settings.N_HORIZON+1, ocp.dims.nx);
u_traj = zeros(settings.N_HORIZON, ocp.dims.nu);

for i=0:settings.N_HORIZON
    x_traj(i+1, :) = ocp_solver.get('x', i);
end
for i=0:settings.N_HORIZON-1
    u_traj(i+1, :) = ocp_solver.get('u', i);
end

% plot trajectories in subplots
t_grid = ocp.solver_options.shooting_nodes;
figure;

subplot(3, 1, 1);
plot(t_grid, x_traj(:, 1), '-');
ylabel('position');
xlim([t_grid(1), t_grid(end)]);

subplot(3, 1, 2);
plot(t_grid, x_traj(:, 2), '-');
ylabel('velocity');
xlim([t_grid(1), t_grid(end)]);

subplot(3, 1, 3);
stairs(t_grid, [u_traj; u_traj(end)], '-');
xlim([t_grid(1), t_grid(end)]);
ylabel('control input');
xlabel('t [s]');

if is_octave()
    waitforbuttonpress;
end

acados_inf = get_acados_infty();

N_HORIZON = 8;    % number of shooting intervals
UMAX = 0.45;

% run
anderson_settings = [0.0, acados_inf, 1e2, 1e1, 1e0];
variant = 'EXACT';
with_abs_cost = true;
tol = 1e-1;

x0 = [0.0; pi; 0.0; 0.0];
Tf = 0.350;
dt_0 = 0.025;

if strcmp(variant, 'GAUSS_NEWTON')
    hessian_approx = 'GAUSS_NEWTON';
    regularize_method = 'NO_REGULARIZE';
    base_label = 'GN';
elseif strcmp(variant, 'EXACT')
    hessian_approx = 'EXACT';
    regularize_method = 'PROJECT';
    base_label = 'project exact Hessian';
else
    error('Unknown variant: %s', variant);
end

solver = setup_ocp_solver(x0, UMAX, dt_0, N_HORIZON, Tf, tol, with_abs_cost, hessian_approx, regularize_method, anderson_settings(1));

t_grid = solver.ocp.solver_options.shooting_nodes;
initial_guess = solver.store_iterate_to_obj();

%
all_kkt_norms = cell(length(anderson_settings), 1);
labels = cell(length(anderson_settings), 1);
sol_list = cell(length(anderson_settings), 1);

for k = 1:numel(anderson_settings)
    anderson_activation_threshold = anderson_settings(k);
    solver.set('anderson_activation_threshold', anderson_activation_threshold);

    % load initial guess, solve and store solution
    solver.load_iterate_from_obj(initial_guess);
    solver.solve();
    solver.print_statistics();
    solution = solver.store_iterate_to_obj();

    % get residuals and compute KKT norms per iteration
    res_all = solver.get('res_all');
    % compute inf-norm across columns for each iteration
    kkt_norms = vecnorm(res_all, inf, 2);
    fprintf('anderson threshold %g -> iterations: %d\n', anderson_activation_threshold, numel(kkt_norms));

    % store results
    if anderson_activation_threshold <= 0.0
        label = base_label;
    elseif acados_inf == anderson_activation_threshold
        label = ['AA(1)-' base_label];
    else
        label = sprintf('AA(1)-%s (thresh=%g)', base_label, anderson_activation_threshold);
    end

    labels{k} = label;
    sol_list{k} = solution;
    all_kkt_norms{k} = kkt_norms;

end

% plot convergence
figure;
for k = 1:numel(all_kkt_norms)
    semilogy(all_kkt_norms{k}, 'DisplayName', labels{k});
    hold on;
end
xlabel('Iteration');
ylabel('KKT norm');
title('Convergence of Furuta pendulum OCP solver');
legend();

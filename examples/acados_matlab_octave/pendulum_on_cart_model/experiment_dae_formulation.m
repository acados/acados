%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.



clear all

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end
import casadi.*

tic

N = 40;

ncases = 3;
constr_violation = zeros(1, ncases);
constr_vals = zeros(N, ncases);
for i = 1:3
    %% arguments
    compile_interface = 'auto';

    gnsf_detect_struct = 'true';

    % discretization
    h = 0.02;

    nlp_solver = 'sqp';
    %nlp_solver = 'sqp_rti';
    nlp_solver_exact_hessian = 'false';
    regularize_method = 'no_regularize';
    nlp_solver_max_iter = 100;
    tol = 1e-12;
    nlp_solver_tol_stat = tol;
    nlp_solver_tol_eq   = tol;
    nlp_solver_tol_ineq = tol;
    nlp_solver_tol_comp = tol;
    nlp_solver_ext_qp_res = 1;
    qp_solver = 'partial_condensing_hpipm';
%     qp_solver = 'full_condensing_hpipm';
%     qp_solver = 'full_condensing_qpoases';
    qp_solver_cond_N = 5;
    qp_solver_cond_ric_alg = 0;
    qp_solver_ric_alg = 0;
    qp_solver_warm_start = 1;
    sim_method = 'irk';
    sim_method_num_stages = 1;
    sim_method_num_steps = 1;
    sim_method_exact_z_output = 0;

    % cost_type = 'linear_ls';
    cost_type = 'ext_cost';
    model_name = ['ocp_pendulum_' num2str(i)];


    %% create model entries
    switch i
        case 1
            model = get_pendulum_on_cart_model();
            theta = model.x(2);
            omega = model.x(4);
            model.con_h_expr = cos(theta)*sin(theta)*omega.^2;
            lh = -40;
            uh = 40;
        case {2,3}
            model = get_pendulum_on_cart_model('dae');
            model.con_h_expr = model.z;
            lh = -40;
            uh = 40;
            if i == 2
                sim_method_exact_z_output = 1;
            end
    end

    % dims
    T = N*h; % horizon length time
    nx = length(model.x);
    nu = length(model.u);

    % constraints
    x0 = [0; pi; 0; 0];
    nbu = nu;
    Jbu = zeros(nbu, nu); for ii=1:nbu; Jbu(ii,ii)=1.0; end
    lbu = -80*ones(nu, 1);
    ubu =  80*ones(nu, 1);


    %% OCP formulation
    ocp = AcadosOcp();
    ocp.model = model;
    ocp.model.name = model_name;
    W_x = diag([1e3, 1e3, 1e-2, 1e-2]);
    W_u = 1e-2;
    ocp.cost.cost_type = 'EXTERNAL';
    ocp.cost.cost_type_e = 'EXTERNAL';
    ocp.model.cost_expr_ext_cost = 0.5 * model.x' * W_x * model.x + 0.5 * model.u' * W_u * model.u;
    ocp.model.cost_expr_ext_cost_e = 0.5 * model.x' * W_x * model.x;
    ocp.model.con_h_expr = model.con_h_expr;
    ocp.constraints.lh = lh;
    ocp.constraints.uh = uh;
    ocp.constraints.x0 = x0;
    ocp.constraints.idxbu = (0:nu-1)';
    ocp.constraints.lbu = lbu;
    ocp.constraints.ubu = ubu;

    ocp.solver_options.N_horizon = N;
    ocp.solver_options.tf = T;
    ocp.solver_options.nlp_solver_type = upper(nlp_solver);
    if strcmp(nlp_solver_exact_hessian, 'true')
        ocp.solver_options.hessian_approx = 'EXACT';
    else
        ocp.solver_options.hessian_approx = 'GAUSS_NEWTON';
    end
    ocp.solver_options.regularize_method = upper(regularize_method);
    ocp.solver_options.nlp_solver_ext_qp_res = nlp_solver_ext_qp_res;
    ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter;
    ocp.solver_options.nlp_solver_tol_stat = nlp_solver_tol_stat;
    ocp.solver_options.nlp_solver_tol_eq = nlp_solver_tol_eq;
    ocp.solver_options.nlp_solver_tol_ineq = nlp_solver_tol_ineq;
    ocp.solver_options.nlp_solver_tol_comp = nlp_solver_tol_comp;
    ocp.solver_options.qp_solver = upper(qp_solver);
    ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N;
    ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg;
    ocp.solver_options.qp_solver_cond_ric_alg = qp_solver_cond_ric_alg;
    ocp.solver_options.qp_solver_warm_start = qp_solver_warm_start;
    ocp.solver_options.integrator_type = 'IRK';
    ocp.solver_options.sim_method_num_stages = sim_method_num_stages;
    ocp.solver_options.sim_method_num_steps = sim_method_num_steps;
    ocp_solver = AcadosOcpSolver(ocp);

    % set trajectory initialization
    x_traj_init = [linspace(0, 0, N+1); linspace(pi, 0, N+1); linspace(0, 0, N+1); linspace(0, 0, N+1)];
    u_traj_init = zeros(nu, N);

    ocp_solver.set('init_x', x_traj_init);
    ocp_solver.set('init_u', u_traj_init);

    ocp_solver.solve();

    % get solution
    u = ocp_solver.get('u');
    x = ocp_solver.get('x');

    %% evaluation
    status = ocp_solver.get('status');
    sqp_iter = ocp_solver.get('sqp_iter');
    time_tot = ocp_solver.get('time_tot');
    time_lin = ocp_solver.get('time_lin');
    time_reg = ocp_solver.get('time_reg');
    time_qp_sol = ocp_solver.get('time_qp_sol');

    fprintf(['\nstatus = %d, sqp_iter = %d, time_int = %f [ms]'...
        ' (time_lin = %f [ms], time_qp_sol = %f [ms], time_reg = %f [ms])\n'],...
        status, sqp_iter, time_tot*1e3, time_lin*1e3, time_qp_sol*1e3, time_reg*1e3 );

    ocp_solver.print('stat');

    if i == 1
        % save reference
        xref = x;
        uref = u;
    else
        % compare error w.r.t reference
        err_x(i) = norm(x - xref)
        err_u(i) = norm(u - uref)
    end

    % compare z accuracy to respective value obtained from x
    thetas = xref(2,:);
    omegas = xref(4,:);
    z_fromx = cos(thetas).*sin(thetas).*omegas.^2;
    if i > 1
        z = ocp_solver.get('z');
        err_z_zfromx(i) = norm( z - z_fromx(1:end-1) );
    end

    % check constraint violation
    theta = model.x(2);
    omega = model.x(4);
    constr_expr = cos(theta)*sin(theta)*omega.^2;
    if i > 1
        z = model.z;
    else
        z = SX.sym('z');
    end
    constr_fun = Function('constr_fun', {model.x, model.u, z}, ...
        {constr_expr});

    constr_violation(i) = 0;
    for j=1:N
        valh = full( constr_fun(x(:,j), u(:,j), z_fromx(:,j) ) );
        constr_vals(j,i) = valh;
        violation = max([0, -(valh-lh), valh - uh]);
        constr_violation(i) = max(norm(violation), constr_violation(i));
    end
end

toc

fprintf('\nConstraint values\n');
fprintf('\nODE \t\tDAE-exact_z\tDAE-extrapolation\n');
for i = 1:N
    fprintf('%.4e\t%.4e\t%.4e\n', constr_vals(i,1), constr_vals(i,2),...
        constr_vals(i,3))
end
fprintf('\n\nConstraint violations');
fprintf('\nODE \t\tDAE-exact_z\tDAE-extrapolation\n');
fprintf('%.4e\t%.4e\t%.4e\n', constr_violation(1), constr_violation(2),...
    constr_violation(3))

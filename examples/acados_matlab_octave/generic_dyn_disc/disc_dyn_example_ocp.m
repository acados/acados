%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.


clear all; clc;
check_acados_requirements()
addpath('../linear_mass_spring_model/');

T = 10.0; % horizon length time
N = 20;
shooting_nodes = linspace(0,10,N+1);
tol = 1e-10;
compare_tol = 1e-8;

variants = {
    struct('name', 'casadi dynamics/casadi cost', 'casadi_dynamics', true, 'casadi_cost', true)
    struct('name', 'generic dynamics/casadi cost', 'casadi_dynamics', false, 'casadi_cost', true)
    struct('name', 'casadi dynamics/generic cost', 'casadi_dynamics', true, 'casadi_cost', false)
    struct('name', 'generic dynamics/generic cost', 'casadi_dynamics', false, 'casadi_cost', false)
};

results = struct('name', {}, 'x', {}, 'u', {}, 'status', {}, 'sqp_iter', {});

for i = 1:length(variants)
    variant = variants{i};

    model = get_linear_mass_spring_model();
    nx = length(model.x);
    nu = length(model.u);
    sym_x = model.x;
    sym_u = model.u;

    ocp = AcadosOcp();
    ocp.model = model;
    ocp.model.name = sprintf('lin_mass_%d_%d', variant.casadi_cost, variant.casadi_dynamics);

    ocp.solver_options.tf = T;
    ocp.solver_options.N_horizon = N;
    ocp.solver_options.shooting_nodes = shooting_nodes;
    ocp.solver_options.nlp_solver_type = 'SQP';
    ocp.solver_options.hessian_approx = 'EXACT';
    ocp.solver_options.nlp_solver_ext_qp_res = 0;
    ocp.solver_options.nlp_solver_max_iter = 100;
    ocp.solver_options.nlp_solver_tol_stat = tol;
    ocp.solver_options.nlp_solver_tol_eq = tol;
    ocp.solver_options.nlp_solver_tol_ineq = tol;
    ocp.solver_options.nlp_solver_tol_comp = tol;
    ocp.solver_options.print_level = 2;
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';
    ocp.solver_options.qp_solver_cond_N = 5;
    ocp.solver_options.integrator_type = 'DISCRETE';

    ocp.cost.cost_type_0 = 'EXTERNAL';
    ocp.cost.cost_type = 'EXTERNAL';
    ocp.cost.cost_type_e = 'EXTERNAL';

    yr_u = zeros(nu, 1);
    yr_x = zeros(nx, 1);
    dWu = 2 * ones(nu, 1);
    dWx = ones(nx, 1);

    cost_expr_ext_cost_0 = 0.5 * sym_u' * (dWu .* sym_u);
    cost_expr_ext_cost = 0.5 * ([sym_u; sym_x] - [yr_u; yr_x])' * ([dWu; dWx] .* ([sym_u; sym_x] - [yr_u; yr_x]));
    cost_expr_ext_cost_e = 0.5 * (sym_x - yr_x)' * (dWx .* (sym_x - yr_x));

    if variant.casadi_cost
        ocp.model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
        ocp.model.cost_expr_ext_cost = cost_expr_ext_cost;
        ocp.model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
    else
        ocp.cost.cost_ext_fun_type_0 = 'generic';
        ocp.cost.cost_generic_source_0 = 'generic_ext_cost.c';
        ocp.cost.cost_generic_fun_0 = 'ext_cost_0_fun';
        ocp.cost.cost_generic_fun_jac_0 = 'ext_cost_0_fun_jac';
        ocp.cost.cost_generic_fun_jac_hess_0 = 'ext_cost_0_fun_jac_hess';

        ocp.cost.cost_ext_fun_type = 'generic';
        ocp.cost.cost_generic_source = 'generic_ext_cost.c';
        ocp.cost.cost_generic_fun = 'ext_cost_fun';
        ocp.cost.cost_generic_fun_jac = 'ext_cost_fun_jac';
        ocp.cost.cost_generic_fun_jac_hess = 'ext_cost_fun_jac_hess';

        ocp.cost.cost_ext_fun_type_e = 'generic';
        ocp.cost.cost_generic_source_e = 'generic_ext_cost.c';
        ocp.cost.cost_generic_fun_e = 'ext_cost_e_fun';
        ocp.cost.cost_generic_fun_jac_e = 'ext_cost_e_fun_jac';
        ocp.cost.cost_generic_fun_jac_hess_e = 'ext_cost_e_fun_jac_hess';
    end

    if variant.casadi_dynamics
        ocp.model.dyn_ext_fun_type = 'casadi';
    else
        ocp.model.disc_dyn_expr = model.disc_dyn_expr;
        ocp.model.dyn_ext_fun_type = 'generic';
        ocp.model.dyn_generic_source = 'generic_disc_dyn.c';
        ocp.model.dyn_disc_fun = 'disc_dyn_fun';
        ocp.model.dyn_disc_fun_jac = 'disc_dyn_fun_jac';
        ocp.model.dyn_disc_fun_jac_hess = 'disc_dyn_fun_jac_hess';
    end

    constr_expr_h = [sym_u; sym_x];
    constr_expr_h_e = sym_x;
    model.con_h_expr = constr_expr_h;
    model.con_h_expr_e = constr_expr_h_e;
    ocp.model.con_h_expr = constr_expr_h;
    ocp.model.con_h_expr_e = constr_expr_h_e;

    x0 = zeros(nx, 1);
    x0(1) = 2.5;
    x0(2) = 2.5;
    ocp.constraints.x0 = x0;
    ocp.constraints.lh = -[0.5 * ones(nu, 1); 4.0 * ones(nx, 1)];
    ocp.constraints.uh = [0.5 * ones(nu, 1); 4.0 * ones(nx, 1)];
    ocp.constraints.lh_e = -4.0 * ones(nx, 1);
    ocp.constraints.uh_e = 4.0 * ones(nx, 1);

    ocp_solver = AcadosOcpSolver(ocp);
    ocp_solver.set('constr_x0', x0);
    ocp_solver.set('init_x', zeros(nx, N+1));
    ocp_solver.set('init_u', zeros(nu, N));

    ocp_solver.solve();
    status = ocp_solver.get('status');
    utraj = ocp_solver.get('u');
    xtraj = ocp_solver.get('x');
    sqp_iter = ocp_solver.get('sqp_iter');

    if status ~= 0
        error('Variant "%s" returned status %d', variant.name, status);
    elseif sqp_iter > 2
        error('Variant "%s" required %d SQP iterations', variant.name, sqp_iter);
    end

    results(i).name = variant.name;
    results(i).x = xtraj;
    results(i).u = utraj;
    results(i).status = status;
    results(i).sqp_iter = sqp_iter;

    fprintf('\n%s: success, sqp_iter = %d\n', variant.name, sqp_iter);
    clear ocp_solver
end

reference = results(1);
fprintf('\nComparison against reference: Casadi/Casadi\n');
for i = 2:length(results)
    err_x = max(abs(reference.x(:) - results(i).x(:)));
    err_u = max(abs(reference.u(:) - results(i).u(:)));
    fprintf('%s: max |x| error = %.3e, max |u| error = %.3e\n', results(i).name, err_x, err_u);

    if err_x > compare_tol || err_u > compare_tol
        error('Variant "%s" differs from the Casadi/Casadi reference', results(i).name);
    end
end

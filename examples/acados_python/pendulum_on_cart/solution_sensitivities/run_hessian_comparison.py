from test_solution_sens_and_exact_hess import *


def run_hessian_comparison(linearized_dynamics=False, discrete=False):
    casadi_solver, lag_hess_fun, lbg, ubg = create_casadi_solver(
        linearized_dynamics=linearized_dynamics,
        discrete=discrete,
    )

    acados_ocp_solver_gn = create_solver(
        hessian_approx='GAUSS_NEWTON',
        linearized_dynamics=linearized_dynamics,
        discrete=discrete,
    )
    acados_ocp_solver_exact = create_solver(
        hessian_approx='EXACT',
        linearized_dynamics=linearized_dynamics,
        discrete=discrete,
    )

    x0 = X0.copy()
    x0[1] = 0.1

    nlp_sol = casadi_solver(p=x0, lbg=lbg, ubg=ubg)
    casadi_hess_l = lag_hess_fun(
        x=nlp_sol['x'],
        p=x0,
        lam_f=1.0,
        lam_g=nlp_sol['lam_g'],
    )['triu_hess_gamma_x_x']
    casadi_hess = ca.triu2symm(ca.triu(casadi_hess_l)).full()
    acados_ocp_solver_gn.solve_for_x0(x0)
    iterate = acados_ocp_solver_gn.get_flat_iterate()
    acados_ocp_solver_exact.set_iterate(iterate)
    acados_ocp_solver_exact.solve_for_x0(
        x0,
        fail_on_nonzero_status=False,
        print_stats_on_failure=False,
    )

    _ = compare_hessian(casadi_hess, acados_ocp_solver_exact)


if __name__ == "__main__":
    run_hessian_comparison(linearized_dynamics=False, discrete=True)
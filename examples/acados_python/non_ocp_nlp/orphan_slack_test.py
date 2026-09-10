import casadi as ca
from acados_template import AcadosOcpSolver, AcadosCasadiOcpQpSolver, AcadosOcp, ACADOS_INFTY, latexify_plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def create_solver(solver_name, variant):

    nv = 2
    ns = 2
    zu = 1 * np.ones(ns)
    Zu = 1 * np.ones(ns)
    zl = np.ones(ns)
    Zl = np.ones(ns)
    q_grad = 10 * np.ones((nv,))
    Q_mat = np.diag([1e1, 1e0])

    ocp = AcadosOcp()
    ocp.solver_options.N_horizon = 0
    v = ca.SX.sym('v', nv)
    ocp.model.x = v
    ocp.model.name = "dense"
    ocp.model.cost_expr_ext_cost_e = v.T @ Q_mat @ v - ca.DM(q_grad).T @ v
    ocp.cost.cost_type_e = "EXTERNAL"

    nb = 2
    ocp.constraints.idxbx_e = np.arange(nb)
    ocp.constraints.lbx_e = -np.ones((nb,))
    ocp.constraints.ubx_e = np.ones((nb,))

    # add slacks
    idxs_rev = np.arange(2)

    ocp.constraints.idxs_rev_e = idxs_rev
    ocp.cost.zl_e = zl
    ocp.cost.Zl_e = Zl
    ocp.cost.zu_e = zu
    ocp.cost.Zu_e = Zu

    if variant.startswith("soft_constrained_masked_irrelevant"):
        ocp.constraints.ubx_e = np.array([ACADOS_INFTY, 1.0])
        ocp.cost.zu_e = zu
        ocp.cost.Zu_e = Zu
    if variant == "soft_constrained_masked_irrelevant_mask_slacks":
        zu[0] = 0
        ocp.cost.zu_e = zu
        ocp.constraints.us_e = -np.array([ACADOS_INFTY, 0.0])

    ocp.solver_options.qp_solver = solver_name
    solver = AcadosOcpSolver(ocp, verbose=False)

    return solver

def solve_with_tol(variant: str, solver: AcadosOcpSolver, tol_comp, tol_others):

    solver.reset()

    solver.options_set('qp_tol_stat', tol_others)
    solver.options_set('qp_tol_eq', tol_others)
    solver.options_set('qp_tol_ineq', tol_others)
    solver.options_set('qp_tol_comp', tol_comp)

    solver.options_set('tol_stat', tol_others)
    solver.options_set('tol_eq', tol_others)
    solver.options_set('tol_ineq', tol_others)
    solver.options_set('tol_comp', tol_comp)

    solver.solve()
    sol = solver.get_iterate()

    # solver.print_statistics()
    # print(f"solving {variant} with: tol_comp {tol_comp:.2e} tol_others {tol_others:.2e}, got status {solver._status}")
    
    qp_iter = sum(solver.get_stats('qp_iter'))

    # extract problem info
    ocp = solver.ocp
    nb = ocp.dims.nbx_e
    ns = ocp.dims.ns_e
    # zu = ocp.cost.zu_e
    # Zu = ocp.cost.Zu_e

    lam_l = sol.lam[0][:nb]
    lam_u = sol.lam[0][nb:2*nb]
    lam_sl = sol.lam[0][2*nb:2*nb+ns]
    lam_su = sol.lam[0][2*nb+ns:2*nb+2*ns]

    su = sol.su[0]
    sl = sol.sl[0]

    # print(f"solution for {variant}: {sol}")
    # print(f"slacks sl: {sol.sl[0]}")
    # print(f"slacks su: {sol.su[0]}")
    # print(f"{lam_l=} {lam_u=}, {lam_sl=} {lam_su=}")

    np.testing.assert_allclose(su[0], 0.0, atol=1e-16, err_msg='orphan slack should be very close to zero')
    np.testing.assert_allclose(lam_u[0], 0.0, atol=1e-16, err_msg='multiplier of masked constraint should be very close to zero')
    np.testing.assert_allclose(lam_su[0], 0.0, atol=1e-16, err_msg='multiplier of masked slack bound should be very close to zero')
    np.testing.assert_array_less(-su[1], 0.0, err_msg='su[1] slack should be strictly positive')

    np.testing.assert_allclose(sl, 0.0, atol=max(tol_comp, tol_others), err_msg='sl should be very close to zero')
    np.testing.assert_allclose(lam_l, 0.0, atol=max(tol_comp, tol_others), err_msg='lam_l should be very close to zero')
    np.testing.assert_array_less(-lam_sl, -0.1, err_msg='lam_sl should be strictly positive')

    return sol, qp_iter


def solve_qp(solver_name: str = 'HPIPM',
         variant='unconstrained',
         tol_comp=1e-6,
         tol_other=1e-2):

    solver: AcadosOcpSolver = create_solver(solver_name, variant)
    sol, qp_iter = solve_with_tol(variant, solver, tol_comp, tol_other)    

    return sol, qp_iter


def slack_su_tol_experiment(variant="soft_constrained_masked_irrelevant", with_plots=False):
    solver_name = "PARTIAL_CONDENSING_HPIPM"
    x_values = np.logspace(-1, -7, 10)
    tol_other_values = [1e-3]

    slack_values_by_tol_other = []
    qp_iters_by_tol_other = []

    solver: AcadosOcpSolver = create_solver(solver_name, variant)

    for tol_other in tol_other_values:
        slack_values = []
        qp_iters = []

        for tol_comp in x_values:
            sol, qp_iter = solve_with_tol(variant, solver, tol_comp, tol_other)

            slack_values.append(abs(sol.su[0][0]))
            qp_iters.append(qp_iter)

        slack_values_by_tol_other.append(np.asarray(slack_values))
        qp_iters_by_tol_other.append(np.asarray(qp_iters))

    if with_plots:
        latexify_plot()
        plt.figure(figsize=(6, 4))
        for tol_other, slack_values in zip(tol_other_values, slack_values_by_tol_other):
            plt.loglog(x_values, slack_values, marker='o', label=r'tol$_{\mathrm{other}}=$' + f'${tol_other:.0e}$')
        plt.xlabel('tol_comp')
        plt.ylabel('$|s_{u,0}|$')
        if 'mask_slacks' in variant:
            title = 'with masked slack bounds, only $L_2$ penalty'
        else:
            title = 'with bounded slacks, $L_1$ and $L_2$ penalty'
        plt.title(title)
        plt.legend()
        plt.grid(True, which='both', ls='--', alpha=0.4)
        plt.tight_layout()
        plt.savefig(f'slack_su_tol_experiment_{variant}.png', dpi=200)

        plt.figure(figsize=(6, 4))
        for tol_other, qp_iters in zip(tol_other_values, qp_iters_by_tol_other):
            plt.semilogx(x_values, qp_iters, marker='o', label=r'tol$_{\mathrm{other}}=$' + f'${tol_other:.0e}$')
        plt.xlabel('tol_comp')
        plt.ylabel('qp iterations')
        plt.title('QP iterations vs comp tolerance')
        plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.legend()
        plt.grid(True, which='both', ls='--', alpha=0.4)
        plt.tight_layout()
        plt.savefig('slack_su_tol_experiment_qp_iters.png', dpi=200)


if __name__ == "__main__":
    with_plots = True
    slack_su_tol_experiment("soft_constrained_masked_irrelevant_mask_slacks", with_plots=with_plots)
    slack_su_tol_experiment("soft_constrained_masked_irrelevant", with_plots=with_plots)
    plt.show()

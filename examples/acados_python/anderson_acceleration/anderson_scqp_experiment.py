from dataclasses import dataclass
import casadi as ca
import numpy as np
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, AcadosOcpFlattenedIterate, AcadosCasadiOcpSolver, ACADOS_INFTY, latexify_plot

from typing import Optional

from scqp_test_problem import build_acados_test_problem
import matplotlib.pyplot as plt

def plot_convergence(list_data: list,
                    list_labels: list,
                    xlim: Optional[int] = None,
                    ylim: tuple=None,
                    save_name: str = None):
    latexify_plot()

    assert len(list_data) == len(list_labels), f"Lists of data and labels do not have the same length, got {len(list_data)} and {len(list_labels)}"

    plt.figure(figsize=(4.5, 3.0))
    plt.clf()
    for i in range(len(list_data)):
        iters = np.arange(0, len(list_data[i]))
        data = np.array(list_data[i]).squeeze()
        plt.semilogy(iters, data, label=list_labels[i])
    plt.legend(loc='best')
    plt.xlabel("iteration number")
    plt.ylabel("KKT residual norm")
    if ylim is not None:
        plt.ylim(ylim)
    if xlim is not None:
        plt.xlim(0, xlim)
    else:
        plt.xlim(0, max([len(data) for data in list_data]))
    if save_name is not None:
        plt.savefig(save_name+".pdf", dpi=300, bbox_inches='tight', pad_inches=0.01)
    plt.tight_layout()
    plt.grid()
    plt.show()

@dataclass
class ExperimentAcadosSettings:
    method: str
    with_anderson_acceleration: bool
    globalization: str = "FIXED_STEP"
    max_iter: int = 200

    def get_label(self):
        label = self.method
        if self.with_anderson_acceleration:
            label = 'AA(1)-' + label
        if self.globalization != "FIXED_STEP":
            label = label + '-' + self.globalization
        return label

@dataclass
class ExperimentResults:
    kkt_norms: np.ndarray
    sol: AcadosOcpFlattenedIterate
    xtraj: np.ndarray
    utraj: np.ndarray


def solve_with_acados(settings: ExperimentAcadosSettings,
                      initial_guess: AcadosOcpFlattenedIterate = None):
    ocp = build_acados_test_problem(mode=settings.method,
                                with_anderson_acceleration=settings.with_anderson_acceleration,
                                globalization=settings.globalization,
                                max_iter=settings.max_iter,
                                )
    N_horizon = ocp.solver_options.N_horizon

    acados_solver = AcadosOcpSolver(ocp)
    if initial_guess is not None:
        acados_solver.load_iterate_from_flat_obj(initial_guess)
    else:
        pass

    # solve
    status = acados_solver.solve()
    acados_solver.print_statistics()
    res_all = acados_solver.get_stats('res_all')
    kkt_norms = np.linalg.norm(res_all, axis=1)
    sol = acados_solver.store_iterate_to_flat_obj()
    acados_solver.dump_last_qp_to_json("qp.json", overwrite=True)
    qp_diag = acados_solver.qp_diagnostics()
    print("qp_diag: ", qp_diag)

    # get solution
    xtraj = np.zeros((N_horizon+1, ocp.dims.nx))
    utraj = np.zeros((N_horizon, ocp.dims.nu))
    for i in range(N_horizon):
        xtraj[i,:] = acados_solver.get(i, "x")
        utraj[i,:] = acados_solver.get(i, "u")
    xtraj[N_horizon,:] = acados_solver.get(N_horizon, "x")
    cost = acados_solver.get_cost()
    print("cost: ", cost)

    results = ExperimentResults(kkt_norms=kkt_norms, sol=sol, xtraj=xtraj, utraj=utraj)

    del acados_solver
    return results


def solve_acados_formulation_with_ipopt(initial_guess: AcadosOcpFlattenedIterate = None):
    ocp = build_acados_test_problem(mode="EXACT")

    N_horizon = ocp.solver_options.N_horizon
    casadi_nlp, bounds = AcadosCasadiOcpSolver.create_casadi_nlp_formulation(ocp)

    casadi_solver = ca.nlpsol("nlp_solver", 'ipopt', casadi_nlp)
                            #   {'ipopt': {'print_level': 1}, 'print_time': False})
    if initial_guess is not None:
        x0 = np.concatenate((initial_guess.x, initial_guess.u))
        nlp_sol = casadi_solver(x0=x0, lbx=bounds['lbx'], ubx=bounds['ubx'], lbg=bounds['lbg'], ubg=bounds['ubg'])
    else:
        nlp_sol = casadi_solver(lbx=bounds['lbx'], ubx=bounds['ubx'], lbg=bounds['lbg'], ubg=bounds['ubg'])


    sol_w = nlp_sol['x']
    xtraj = sol_w[:ocp.dims.nx*(N_horizon+1)].reshape((ocp.dims.nx, N_horizon+1)).full()
    utraj = sol_w[ocp.dims.nx*(N_horizon+1):].reshape((ocp.dims.nu, N_horizon)).full()

    xtraj = xtraj.T
    utraj = utraj.T
    results = ExperimentResults(kkt_norms=None, sol=None, xtraj=xtraj, utraj=utraj)
    return results


def main_acados():
    ref_settings = ExperimentAcadosSettings(method='SCQP', with_anderson_acceleration=False, globalization='MERIT_BACKTRACKING')
    ref_res = solve_with_acados(ref_settings)

    # set up initial guess
    sol = ref_res.sol
    perturb_scale = 1e-6
    acados_guess = AcadosOcpFlattenedIterate(
        x = sol.x + perturb_scale * np.random.randn(sol.x.shape[0]),
        u = sol.u + perturb_scale * np.random.randn(sol.u.shape[0]),
        z = sol.z + 0 * np.random.randn(sol.z.shape[0]),
        sl = sol.sl + 0 * np.random.randn(sol.sl.shape[0]),
        su = sol.su + 0 * np.random.randn(sol.su.shape[0]),
        pi = sol.pi, # perturb_scale * np.random.randn(sol.pi.shape[0]),
        lam = np.abs(sol.lam + 1 * np.random.randn(sol.lam.shape[0])),
    )

    ipopt_res = solve_acados_formulation_with_ipopt(initial_guess=sol)

    # setup acados guess to match ipopt solution
    acados_guess.x = ipopt_res.xtraj.flatten()
    acados_guess.u = ipopt_res.utraj.flatten()
    ref_res_2 = solve_with_acados(ref_settings, initial_guess=acados_guess)

    # compare with ipopt
    print("compare: IPOPT vs acados")
    diff_x = np.linalg.norm(ref_res.xtraj - ipopt_res.xtraj)
    diff_u = np.linalg.norm(ref_res.utraj - ipopt_res.utraj)
    print("diff x: ", diff_x)
    print("diff u: ", diff_u)
    # print("diff x: ", ref_res.xtraj - ipopt_res.xtraj)
    # print("diff u: ", ref_res.utraj - ipopt_res.utraj)
    print("compare: IPOPT vs acados started at IPOPT solution")
    diff_x = np.linalg.norm(ref_res_2.xtraj - ipopt_res.xtraj)
    diff_u = np.linalg.norm(ref_res_2.utraj - ipopt_res.utraj)
    print("diff x: ", diff_x)
    print("diff u: ", diff_u)

    # disturb again
    sol = ref_res_2.sol
    perturb_scale = 1e-3
    acados_guess = AcadosOcpFlattenedIterate(
        x = sol.x + perturb_scale * np.random.randn(sol.x.shape[0]),
        u = sol.u + perturb_scale * np.random.randn(sol.u.shape[0]),
        z = sol.z + 0 * np.random.randn(sol.z.shape[0]),
        sl = sol.sl + 0 * np.random.randn(sol.sl.shape[0]),
        su = sol.su + 0 * np.random.randn(sol.su.shape[0]),
        pi = sol.pi, # perturb_scale * np.random.randn(sol.pi.shape[0]),
        lam = np.abs(sol.lam + 1 * np.random.randn(sol.lam.shape[0])),
    )

    # plot solution
    # plot_pendulum(np.linspace(0, ocp_problem.dt*ocp_problem.N, ocp_problem.N+1),
    #               ref_res.utraj, ref_res.xtraj, show=False)
    # plot_pendulum(np.linspace(0, ocp_problem.dt*ocp_problem.N, ocp_problem.N+1),
    #               ipopt_res.utraj, ipopt_res.xtraj, show=False)
    # plot_pendulum(np.linspace(0, ocp_problem.dt*ocp_problem.N, ocp_problem.N+1),
    #             ref_res_2.utraj, ref_res_2.xtraj, show=False)
    # plt.show()

    # Experiment
    settings = [
        ExperimentAcadosSettings(method='EXACT', with_anderson_acceleration=False, max_iter=100),
        ExperimentAcadosSettings(method='GN', with_anderson_acceleration=False, max_iter=60),
        ExperimentAcadosSettings(method='SCQP', with_anderson_acceleration=False, max_iter=60),
        ExperimentAcadosSettings(method='SCQP', with_anderson_acceleration=True),
        # ExperimentAcadosSettings(method='GN', with_anderson_acceleration=True),
    ]
    # Evaluation
    labels = [s.get_label() for s in settings]
    results = []
    for i, setting in enumerate(settings):
        res = solve_with_acados(setting, initial_guess=acados_guess)
        results.append(res)
    plot_convergence([r.kkt_norms for r in results], labels)


if __name__ == "__main__":
    main_acados()

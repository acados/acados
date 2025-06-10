#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#


from dataclasses import dataclass
import numpy as np
from acados_template import AcadosOcpSolver, AcadosOcpFlattenedIterate, AcadosCasadiOcpSolver, plot_convergence

from scqp_test_problem import build_acados_test_problem, plot_pendulum

import matplotlib.pyplot as plt

np.random.seed(0)

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


def solve_with_acados(settings: ExperimentAcadosSettings,
                      initial_guess: AcadosOcpFlattenedIterate = None):
    ocp = build_acados_test_problem(mode=settings.method,
                                with_anderson_acceleration=settings.with_anderson_acceleration,
                                globalization=settings.globalization,
                                max_iter=settings.max_iter,
                                )
    acados_solver = AcadosOcpSolver(ocp, verbose=False)
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
    # qp_diag = acados_solver.qp_diagnostics()
    # print("qp_diag: ", qp_diag)

    cost = acados_solver.get_cost()
    print("cost: ", cost)

    results = ExperimentResults(kkt_norms=kkt_norms, sol=sol)

    del acados_solver
    return results


def solve_acados_formulation_with_ipopt(initial_guess: AcadosOcpFlattenedIterate = None,
                                        mode="EXACT"):
    ocp = build_acados_test_problem(mode=mode)

    solver = AcadosCasadiOcpSolver(ocp, use_acados_hessian=True)

    if initial_guess is not None:
        solver.load_iterate_from_flat_obj(initial_guess)
    solver.solve()

    sol = solver.store_iterate_to_flat_obj()

    results = ExperimentResults(kkt_norms=None, sol=sol)
    return results

def raise_test_failure_message(msg: str):
    # print(f"ERROR: {msg}")
    raise Exception(msg)

def main(plot_sol=False):
    ref_settings = ExperimentAcadosSettings(method='SCQP', with_anderson_acceleration=False, globalization='MERIT_BACKTRACKING')
    ref_res = solve_with_acados(ref_settings)
    sol = ref_res.sol

    # ipopt_res = solve_acados_formulation_with_ipopt(initial_guess=sol, mode='EXACT')
    # NOTE: the above converges to a worse local optimum.
    ipopt_res = solve_acados_formulation_with_ipopt(initial_guess=sol, mode='SCQP')

    if plot_sol:
        ocp = build_acados_test_problem()
        ocp.make_consistent()
        xtraj = ipopt_res.sol.x.reshape((ocp.dims.N+1, ocp.dims.nx))
        utraj = ipopt_res.sol.u.reshape((ocp.dims.N, ocp.dims.nu))
        plot_pendulum(ocp.solver_options.shooting_nodes,
                      utraj, xtraj, plt_show=False)
        plt.show()

    # start acados at ipopt solution
    ref_res_2 = solve_with_acados(ref_settings, initial_guess=ipopt_res.sol)

    # compare with ipopt
    print("compare: IPOPT vs acados")
    diff_x = np.linalg.norm(ref_res.sol.x - ipopt_res.sol.x)
    diff_u = np.linalg.norm(ref_res.sol.u - ipopt_res.sol.u)
    print("diff x: ", diff_x)
    print("diff u: ", diff_u)
    # print("diff x: ", ref_res.sol.x - ipopt_res.sol.x)
    # print("diff u: ", ref_res.sol.u - ipopt_res.sol.u)
    print("compare: IPOPT vs acados started at IPOPT solution")
    diff_x = np.linalg.norm(ref_res_2.sol.x - ipopt_res.sol.x)
    diff_u = np.linalg.norm(ref_res_2.sol.u - ipopt_res.sol.u)
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
    results: list[ExperimentResults] = []
    for i, setting in enumerate(settings):
        res = solve_with_acados(setting, initial_guess=acados_guess)
        results.append(res)
    plot_convergence([r.kkt_norms for r in results], labels)

    # asserts to check behavior
    res_anderson = results[-1]
    res_exact = results[0]
    res_scqp = results[-2]

    def assert_convergence_iterations(res, min_iter: int, max_iter: int, method: str):
        n_iter = len(res.kkt_norms)
        if n_iter > max_iter or n_iter < min_iter:
            raise_test_failure_message(f"Number of iterations {n_iter} not in expected range [{min_iter}, {max_iter}] for {method}")

    # assert all solutions are the same
    ref_sol = results[0].sol
    for i, (res, setting) in enumerate(zip(results[1:], settings[1:])):
        if setting.method == "GN":
            print(f"Skipping GN comparison for {setting.get_label()} expect non converged result.")
            continue
        diff_x = np.linalg.norm(ref_sol.x - res.sol.x)
        diff_u = np.linalg.norm(ref_sol.u - res.sol.u)
        print(f"diff x for {labels[i+1]}: ", diff_x)
        print(f"diff u for {labels[i+1]}: ", diff_u)
        if diff_x > 1e-6 or diff_u > 1e-6:
            raise_test_failure_message(f"Solution mismatch for {labels[i+1]}: x diff {diff_x}, u diff {diff_u}")
        else:
            print(f"Solution for {labels[i+1]} matches reference solution.")

    assert_convergence_iterations(res_anderson, 8, 12, "Anderson")
    assert_convergence_iterations(res_exact, 4, 6, "Exact")
    assert_convergence_iterations(res_scqp, 25, 40, "SCQP")

if __name__ == "__main__":
    main()

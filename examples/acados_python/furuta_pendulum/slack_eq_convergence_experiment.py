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

from furuta_common import setup_ocp_solver
from utils import plot_furuta_pendulum
import numpy as np

from acados_template import AcadosOcpFlattenedIterate, plot_convergence, plot_contraction_rates
from typing import Tuple
import matplotlib.pyplot as plt
N_HORIZON = 8   # number of shooting intervals
UMAX = .45


def test_solver(with_anderson_acceleration: bool,
                variant='GAUSS_NEWTON') -> Tuple[AcadosOcpFlattenedIterate, np.ndarray]:
    x0 = np.array([0.0, np.pi, 0.0, 0.0])

    Tf = .350       # total prediction time
    dt_0 = 0.025    # sampling time = length of first shooting interval

    if variant == 'GAUSS_NEWTON':
        hessian_approx = 'GAUSS_NEWTON'
        regularize_method = 'NO_REGULARIZE'
    elif variant == 'EXACT':
        hessian_approx = 'EXACT'
        regularize_method = 'PROJECT'
    solver = setup_ocp_solver(x0, UMAX, dt_0, N_HORIZON, Tf, with_anderson_acceleration=with_anderson_acceleration, nlp_solver_max_iter = 500, tol = 1e-8, with_abs_cost=True, hessian_approx=hessian_approx, regularize_method=regularize_method)

    t_grid = solver.acados_ocp.solver_options.shooting_nodes

    status = solver.solve()
    solver.print_statistics()
    solution = solver.store_iterate_to_flat_obj()

    res_all = solver.get_stats('res_all')
    kkt_norms = np.linalg.norm(res_all, axis=1)

    return solution, kkt_norms, t_grid

def raise_test_failure_message(msg: str):
    # print(f"ERROR: {msg}")
    raise Exception(msg)

def main(variant: str='GAUSS_NEWTON', plot_trajectory: bool=False, store_plots: bool=False):
    # test with anderson acceleration
    kkt_norm_list = []
    contraction_rates_list = []
    sol_list = []
    labels = []
    if variant == 'GAUSS_NEWTON':
        base_label = "GN"
    elif variant == 'EXACT':
        base_label = "project exact Hessian"
    else:
        raise ValueError(f"Unknown variant: {variant}")

    for with_anderson_acceleration in [True, False]:
        sol, kkt_norms, t_grid = test_solver(with_anderson_acceleration=with_anderson_acceleration, variant=variant)
        # compute contraction rates
        contraction_rates = kkt_norms[1:-1]/kkt_norms[0:-2]
        # append results
        kkt_norm_list.append(kkt_norms)
        contraction_rates_list.append(contraction_rates)
        sol_list.append(sol)
        label = "AA(1)-" + base_label if with_anderson_acceleration else base_label
        labels.append(label)
        # checks
        n_iter = len(kkt_norms)
        if with_anderson_acceleration:
            if not n_iter < 30:
                raise_test_failure_message(f"Expected less than 30 iterations with Anderson acceleration, got {n_iter}")
        else:
            if not n_iter > 60:
                raise_test_failure_message(f"Expected more than 60 iterations without Anderson acceleration, got {n_iter}")

    # checks
    ref_sol = sol_list[0]
    for i, sol in enumerate(sol_list[1:]):
        if not ref_sol.allclose(sol, atol=1e-4):
            print(f"Solution mismatch for {labels[i]}: difference: {(ref_sol-sol).inf_norm()}")

        else:
            print(f"Solution for {labels[i]} matches reference solution.")

    # plot results
    plot_convergence(
        kkt_norm_list,
        labels,
        fig_filename=f"convergence_{variant}_furuta.png" if store_plots else None
    )
    plot_contraction_rates(
        contraction_rates_list,
        labels,
        fig_filename=f"contraction_rates_{variant}_furuta.png" if store_plots else None
    )
    if plot_trajectory:
        plot_furuta_pendulum(t_grid, ref_sol.x.reshape((N_HORIZON + 1, -1)), ref_sol.u.reshape((N_HORIZON, -1)), UMAX, plt_show=False)
        sol1 = sol_list[1]
        plot_furuta_pendulum(t_grid, sol1.x.reshape((N_HORIZON + 1, -1)), sol1.u.reshape((N_HORIZON, -1)), UMAX, plt_show=False)
        plt.show()

if __name__ == "__main__":
    plot_trajectory = False
    store_plots = False
    main("GAUSS_NEWTON", plot_trajectory=plot_trajectory, store_plots=store_plots)
    main("EXACT", plot_trajectory=plot_trajectory, store_plots=store_plots)


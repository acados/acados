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

from acados_template import (
    AcadosOcpSolver,
    AcadosOcpFlattenedIterate,
    plot_convergence,
    plot_contraction_rates,
    ACADOS_INFTY,
)
from typing import Tuple
import matplotlib.pyplot as plt

N_HORIZON = 8  # number of shooting intervals
UMAX = 0.45


def create_solver(variant, tol, with_abs_cost):
    x0 = np.array([0.0, np.pi, 0.0, 0.0])

    Tf = 0.350  # total prediction time
    dt_0 = 0.025  # sampling time = length of first shooting interval

    if variant == "GAUSS_NEWTON":
        hessian_approx = "GAUSS_NEWTON"
        regularize_method = "NO_REGULARIZE"
    elif variant == "EXACT":
        hessian_approx = "EXACT"
        regularize_method = "PROJECT"
    solver = setup_ocp_solver(
        x0,
        UMAX,
        dt_0,
        N_HORIZON,
        Tf,
        nlp_solver_max_iter=500,
        tol=tol,
        with_abs_cost=with_abs_cost,
        hessian_approx=hessian_approx,
        regularize_method=regularize_method,
        with_anderson_acceleration=True
    )

    return solver


def test_solver(
    solver: AcadosOcpSolver,
    initial_guess: AcadosOcpFlattenedIterate,
    anderson_activation_threshold: float,
):
    solver.options_set("anderson_activation_threshold", anderson_activation_threshold)

    solver.load_iterate_from_flat_obj(initial_guess)
    status = solver.solve()
    solver.print_statistics()
    solution = solver.store_iterate_to_flat_obj()

    res_all = solver.get_stats("res_all")
    kkt_norms = np.linalg.norm(res_all, axis=1)
    return solution, kkt_norms


def raise_test_failure_message(msg: str):
    # print(f"ERROR: {msg}")
    raise Exception(msg)


def main(
    anderson_settings: list,
    variant: str = "GAUSS_NEWTON",
    plot_trajectory: bool = False,
    store_plots: bool = False,
    plot_contraction: bool = False,
    ignore_checks: bool = False,
    with_abs_cost: bool = True,
    tol: float = 1e-8,
):
    # test with anderson acceleration
    kkt_norm_list = []
    contraction_rates_list = []
    sol_list = []
    labels = []
    if variant == "GAUSS_NEWTON":
        base_label = "GN"
    elif variant == "EXACT":
        base_label = "project exact Hessian"
    else:
        raise ValueError(f"Unknown variant: {variant}")

    solver = create_solver(variant, tol, with_abs_cost)

    t_grid = solver.acados_ocp.solver_options.shooting_nodes
    initial_guess = solver.store_iterate_to_flat_obj()

    for anderson_activation_threshold in anderson_settings:

        sol, kkt_norms = test_solver(
            solver, initial_guess, anderson_activation_threshold
        )
        # compute contraction rates
        contraction_rates = kkt_norms[1:-1] / kkt_norms[0:-2]
        # append results
        kkt_norm_list.append(kkt_norms)
        contraction_rates_list.append(contraction_rates)
        sol_list.append(sol)
        if anderson_activation_threshold <= 0.0:
            label = base_label
        elif anderson_activation_threshold == ACADOS_INFTY:
            label = "AA(1)-" + base_label
        else:
            label = f"AA(1)-{base_label} (thresh={anderson_activation_threshold})"
        labels.append(label)

        # checks
        n_iter = len(kkt_norms)
        if ignore_checks:
            print("Ignoring test checks.")
            continue
        if tol < 1e-5:
            if with_abs_cost:
                if anderson_activation_threshold == ACADOS_INFTY: # full Anderson
                    if not n_iter < 30:
                        raise_test_failure_message(
                            f"Expected less than 30 iterations with Anderson acceleration, got {n_iter}"
                        )
                elif anderson_activation_threshold <= 0.0: # no Anderson
                    if not n_iter > 60:
                        raise_test_failure_message(
                            f"Expected more than 60 iterations without Anderson acceleration, got {n_iter}"
                        )
            else:
                if anderson_activation_threshold == ACADOS_INFTY: # full Anderson
                    if not n_iter < 27:
                        raise_test_failure_message(
                            f"Expected less than 27 iterations with Anderson acceleration, got {n_iter}"
                        )
                elif anderson_activation_threshold <= 0.0: # no Anderson
                    if not n_iter > 200:
                        raise_test_failure_message(
                            f"Expected more than 200 iterations without Anderson acceleration, got {n_iter}"
                        )

    # checks
    ref_sol = sol_list[0]
    for i, sol in enumerate(sol_list[1:]):
        if not ref_sol.allclose(sol, atol=1e-4):
            print(
                f"Solution mismatch for {labels[i]}: difference: {(ref_sol-sol).inf_norm()}"
            )

        else:
            print(f"Solution for {labels[i]} matches reference solution.")

    # plot results
    plot_convergence(
        kkt_norm_list,
        labels,
        figsize=(7.0, 4.0),
        fig_filename=(
            f"convergence_{'slack' if with_abs_cost else ''}_{variant}_{tol}_furuta.png"
            if store_plots
            else None
        ),
        title=f"Furuta pendulum {'slack' if with_abs_cost else ''} OCP with tolerance {tol}",
        show_plot=False,
    )
    if plot_contraction:
        plot_contraction_rates(
            contraction_rates_list,
            labels,
            fig_filename=(
                f"contraction_rates_{'slack' if with_abs_cost else ''}_{variant}_{tol}_furuta.png"
                if store_plots
                else None
            ),
            show_plot=False,
        )
    if plot_trajectory:
        plot_furuta_pendulum(
            t_grid,
            ref_sol.x.reshape((N_HORIZON + 1, -1)),
            ref_sol.u.reshape((N_HORIZON, -1)),
            UMAX,
            plt_show=False,
        )
        sol1 = sol_list[1]
        plot_furuta_pendulum(
            t_grid,
            sol1.x.reshape((N_HORIZON + 1, -1)),
            sol1.u.reshape((N_HORIZON, -1)),
            UMAX,
            plt_show=False,
        )
        plt.show()


if __name__ == "__main__":
    plot_trajectory = False
    store_plots = True
    anderson_settings = [0.0, ACADOS_INFTY, 1e2, 1e1, 1e0]

    main(anderson_settings, "GAUSS_NEWTON", with_abs_cost=False, plot_trajectory=plot_trajectory, store_plots=store_plots)
    main(anderson_settings, "GAUSS_NEWTON", plot_trajectory=plot_trajectory, store_plots=store_plots)
    main(anderson_settings, "EXACT", plot_trajectory=plot_trajectory, store_plots=store_plots)
    # Below case shows that AA with very loose tolerance can slow down convergence.
    main(anderson_settings, "EXACT", plot_trajectory=plot_trajectory, tol=1e-1, ignore_checks=True, store_plots=store_plots)
    plt.show()

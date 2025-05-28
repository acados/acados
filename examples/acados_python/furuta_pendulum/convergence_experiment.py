# -*- coding: future_fstrings -*-
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
import numpy as np

from acados_template import AcadosOcpFlattenedIterate, latexify_plot
from typing import Tuple
import matplotlib.pyplot as plt

def plot_convergence(list_data: list,
                    list_labels: list,
                    fig_filename: str = None):
    latexify_plot()

    assert len(list_data) == len(list_labels), f"Lists of data and labels do not have the same length, got {len(list_data)} and {len(list_labels)}"

    plt.figure(figsize=(4.5, 3.0))
    for i in range(len(list_data)):
        iters = np.arange(0, len(list_data[i]))
        data = np.array(list_data[i]).squeeze()
        plt.semilogy(iters, data, label=list_labels[i])
    plt.legend(loc='best')
    plt.xlabel("iteration number")
    plt.ylabel("KKT residual norm")
    plt.xlim(0, max([len(data) for data in list_data]))
    if fig_filename is not None:
        plt.savefig(fig_filename, dpi=300, bbox_inches='tight', pad_inches=0.01)
    plt.tight_layout()
    plt.grid()

def plot_contraction_rates(list_data: list,
                          list_labels: list,
                          fig_filename: str = None):
    latexify_plot()
    plt.figure(figsize=(4.5, 3.0))
    for rates, label in zip(list_data, list_labels):
        iters = np.arange(0, len(rates))
        plt.plot(iters, rates, label=label)
    plt.legend(loc='best')
    plt.xlabel("iteration number")
    plt.ylabel("empirical contraction rate")
    plt.xlim(0, max([len(data) for data in list_data]))
    if fig_filename is not None:
        plt.savefig(fig_filename, dpi=300, bbox_inches='tight', pad_inches=0.01)
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.grid()


def test_solver(with_anderson_acceleration: bool) -> Tuple[AcadosOcpFlattenedIterate, np.ndarray]:
    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    umax = .45

    Tf = .350       # total prediction time
    N_horizon = 8   # number of shooting intervals
    dt_0 = 0.025    # sampling time = length of first shooting interval

    solver = setup_ocp_solver(x0, umax, dt_0, N_horizon, Tf, with_anderson_acceleration=with_anderson_acceleration, nlp_solver_max_iter = 500, tol = 1e-8)

    status = solver.solve()
    solver.print_statistics()
    solution = solver.store_iterate_to_flat_obj()

    res_all = solver.get_stats('res_all')
    kkt_norms = np.linalg.norm(res_all, axis=1)

    return solution, kkt_norms

def raise_test_failure_message(msg: str):
    # print(f"ERROR: {msg}")
    raise Exception(msg)

def main():
    # test with anderson acceleration
    kkt_norm_list = []
    contraction_rates_list = []
    sol_list = []
    labels = []
    for with_anderson_acceleration in [True, False]:
        sol, kkt_norms = test_solver(with_anderson_acceleration=with_anderson_acceleration)
        # compute contraction rates
        contraction_rates = kkt_norms[1:-1]/kkt_norms[0:-2]
        # append results
        kkt_norm_list.append(kkt_norms)
        contraction_rates_list.append(contraction_rates)
        sol_list.append(sol)
        labels.append("AA(1)-GN" if with_anderson_acceleration else "GN")
        # checks
        n_iter = len(kkt_norms)
        if with_anderson_acceleration:
            assert n_iter < 27, f"Expected less than 27 iterations with Anderson acceleration, got {n_iter}"
        else:
            assert n_iter > 200, f"Expected more than 200 iterations without Anderson acceleration, got {n_iter}"

    # checks
    ref_sol = sol_list[0]
    for i, sol in enumerate(sol_list[1:]):
        if not ref_sol.allclose(sol, atol=1e-6):
            raise_test_failure_message(f"Solution mismatch for {labels[i]}: {sol} vs {ref_sol}")
        else:
            print(f"Solution for {labels[i]} matches reference solution.")

    # plot results
    plot_convergence(
        list_data=kkt_norm_list,
        list_labels=labels,
        # fig_filename="convergence_furuta_pendulum.png"
    )
    plot_contraction_rates(
        list_data=contraction_rates_list,
        list_labels=labels,
        # fig_filename="contraction_rates_furuta_pendulum.png"
    )
    plt.show()

if __name__ == "__main__":
    main()


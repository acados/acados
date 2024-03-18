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

# authors: Katrin Baumgaertner, Jonathan Frey

from furuta_model import get_furuta_model

import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from acados_template import latexify_plot, AcadosSim, AcadosSimSolver

latexify_plot()

X0 = np.array([0, np.pi/2, 0, 0])


def setup_acados_integrator(model, dt, num_stages=4, num_steps=4, integrator_type="ERK",
                            newton_iter=20, newton_tol=1e-10):
    sim = AcadosSim()
    sim.model = model
    sim.solver_options.T = dt
    sim.solver_options.num_stages = num_stages
    sim.solver_options.num_steps = num_steps
    sim.solver_options.integrator_type = integrator_type
    sim.solver_options.newton_iter = newton_iter
    sim.solver_options.newton_tol = newton_tol
    sim.solver_options.collocation_type = 'GAUSS_RADAU_IIA'
    # sim.solver_options.sens_forw = False
    sim.solver_options.sens_adj = False
    sim.solver_options.sens_algebraic = False
    sim.solver_options.sens_hess = False
    sim.solver_options.sim_method_jac_reuse = True

    acados_integrator = AcadosSimSolver(sim, verbose=False)
    return acados_integrator


def simulate(
        integrator: AcadosSimSolver, x0, u_trajectory
):
    Nsim = u_trajectory.shape[1]
    x_trajectory = np.zeros((integrator.acados_sim.dims.nx, Nsim + 1))
    timings = np.zeros((Nsim + 1,))
    x_trajectory[:, 0] = x0
    for i in range(Nsim):
        x_trajectory[:, i + 1] = integrator.simulate(x_trajectory[:, i], u_trajectory[:, i])
        timings[i] = integrator.get('time_tot')
    return x_trajectory, timings

def main(integrator_type = "IRK"):
    SAVE_FIG = True

    Tsim = 1.0
    dt_plant = 0.2

    model = get_furuta_model()

    Nsim = int(Tsim / dt_plant)
    if not (Tsim / dt_plant).is_integer():
        print("WARNING: Tsim / dt_plant should be an integer!")

    integrator = setup_acados_integrator(
        model,
        dt_plant,
        num_stages=8,
        num_steps=100,
        newton_iter=20,
        newton_tol=1e-16,
        integrator_type="IRK",
    )

    nu = 1
    np.random.seed(0)
    u_trajectory = 0.03 * (np.random.random((nu, Nsim)) - 0.5)
    X_exact, _ = simulate(integrator, X0, u_trajectory)
    del integrator

    # store results for plotting
    X_all = [X_exact]
    labels_all = ["exact"]

    plt.figure()

    if integrator_type == "ERK":
        num_stages_vals = [1, 2, 4]
        num_steps_vals = [1, 2, 5, 10, 50, 100]
    else:
        num_stages_vals = [1, 2, 3, 4]
        num_steps_vals = [1, 2, 5, 10, 50]
    markers = ["o", "v", "s", "D", "^", ">", "<", "1", "2", "3", "4"]
    colors = [plt.cm.tab10(i) for i in range(len(num_stages_vals))]

    ## EXPERIMENT
    for i, num_stages in enumerate(num_stages_vals):
        for j, num_steps in enumerate(num_steps_vals):

            label = f"{integrator_type}_stages_{num_stages}_steps_{num_steps}"
            print(f"Running simulation with {label}")

            model = get_furuta_model()

            integrator = setup_acados_integrator(
                model,
                dt_plant,
                num_stages=num_stages,
                num_steps=num_steps,
                integrator_type=integrator_type,
            )

            X, timings_integrator = simulate(integrator, X0, u_trajectory)
            err_x = np.max(np.abs(X - X_exact))

            # store all results
            X_all.append(X)
            labels_all.append(label)

            plt.plot(
                [np.mean(timings_integrator * 1e3)],
                [err_x],
                color=colors[i],
                marker=markers[j],
                label=label,
            )

            print(f"got: err_x: {err_x}")

            del integrator

    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim([1e-10, 1e1])
    plt.ylabel(r"error $\Vert x - x^{\mathrm{exact}}\Vert_{\infty}$")
    plt.xlabel("mean integration time [ms]")
    plt.title(f"{integrator_type}")

    legend_elements = [
        Line2D([0], [0], color=colors[i], lw=4, label="num stages = " + str(num_stages))
        for i, num_stages in enumerate(num_stages_vals)
    ] + [
        Line2D(
            [0],
            [0],
            marker=markers[j],
            lw=0,
            color="k",
            label="num steps = " + str(num_steps),
        )
        for j, num_steps in enumerate(num_steps_vals)
    ]

    plt.legend(handles=legend_elements)
    if SAVE_FIG:
        fig_filename = f"cstr_integrator_experiment_{integrator_type}.pdf"
        plt.savefig(
            fig_filename, bbox_inches="tight", transparent=True, pad_inches=0.05
        )
        print(f"\nstored figure in {fig_filename}")

    plt.show()

    # print failed runs
    for x, l in zip(X_all, labels_all):
        if np.any(np.isnan(x)):
            print(f"Run {l} failed with NaN.")

    print(f"final state {X_exact[:, -1]}")


if __name__ == "__main__":
    main(integrator_type="IRK")
    main(integrator_type="ERK")

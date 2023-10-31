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

from cstr_model import CSTRParameters, setup_cstr_model
from setup_acados_ocp_solver import MpcCSTRParameters
from setup_acados_integrator import setup_acados_integrator
from setup_acados_ocp_solver import setup_acados_ocp_solver
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from main import simulate
from acados_template import latexify_plot

latexify_plot()

def main():
    SAVE_FIG = True
    integrator_type = "ERK"  # 'IRK'

    Tsim = 25
    dt_plant = 0.25  # [min]

    cstr_params = CSTRParameters()
    mpc_params = MpcCSTRParameters(xs=cstr_params.xs, us=cstr_params.us)
    model = setup_cstr_model(cstr_params)

    Nsim = int(Tsim / dt_plant)
    if not (Tsim / dt_plant).is_integer():
        print("WARNING: Tsim / dt_plant should be an integer!")

    # steady-state
    xs = np.array([[0.878, 324.5, 0.659]]).T
    us = np.array([[300, 0.1]]).T

    # constant ref
    X_ref = np.tile(xs, Nsim + 1).T
    U_ref = np.tile(us, Nsim).T

    # reference jump
    xs2 = np.array([0.7, 337, 0.75])
    us2 = np.array([305, 0.1])

    Njump = int(Nsim / 3)
    X_ref[Njump : 2 * Njump, :] = xs2
    U_ref[Njump : 2 * Njump, :] = us2

    # initial state
    x0 = np.array([0.05, 0.75, 0.5]) * xs.ravel()

    # compute exact solution
    print("\n\nRunning simulation with reference input\n\n")
    mpc_params = MpcCSTRParameters(xs=cstr_params.xs, us=cstr_params.us)
    ocp_solver = setup_acados_ocp_solver(model, mpc_params, cstr_params=cstr_params)

    integrator = setup_acados_integrator(
        model,
        dt_plant,
        cstr_param=cstr_params,
        num_stages=8,
        num_steps=100,
        integrator_type="IRK",
    )
    X_exact, U_exact, _, _ = simulate(ocp_solver, integrator, x0, Nsim, X_ref, U_ref)

    # use just U_ref
    # X_exact, U_exact, _, _ = simulate(None, integrator, x0, Nsim, X_ref, U_ref)
    del integrator

    # store results for plotting
    X_all = [X_exact]
    U_all = [U_exact]
    labels_all = ["exact"]

    plt.figure()

    if integrator_type == "ERK":
        num_stages_vals = [1, 2, 4]
        num_steps_vals = [1, 2, 5, 10, 50, 100, 500]
    else:
        num_stages_vals = [1, 2, 3, 4, 5]
        num_steps_vals = [1, 2, 5, 10, 50, 100]
    markers = ["o", "v", "s", "D", "^", ">", "<", "1", "2", "3", "4"]
    colors = [plt.cm.tab10(i) for i in range(len(num_stages_vals))]

    ## EXPERIMENT
    for i, num_stages in enumerate(num_stages_vals):
        for j, num_steps in enumerate(num_steps_vals):

            label = f"{integrator_type}_stages_{num_stages}_steps_{num_steps}"
            print(f"\n\nRunning simulation with {label}\n\n")

            model = setup_cstr_model(cstr_params)
            # integrator = setup_acados_integrator(model, dt_plant, cstr_param=cstr_params, num_stages=num_stages, num_steps=num_steps, integrator_type=integrator_type, name=label)
            integrator = setup_acados_integrator(
                model,
                dt_plant,
                cstr_param=cstr_params,
                num_stages=num_stages,
                num_steps=num_steps,
                integrator_type=integrator_type,
            )

            # X, U, _, timings_integrator = simulate(None, integrator, x0, Nsim, X_ref, U_ref)
            X, U, _, timings_integrator = simulate(
                None, integrator, x0, Nsim, X_exact, U_exact
            )
            err_x = np.max(np.abs(X - X_exact))
            err_u = np.max(np.abs(U - U_exact))

            # store all results
            X_all.append(X)
            U_all.append(U)
            labels_all.append(label)

            plt.plot(
                [np.mean(timings_integrator * 1e3)],
                [err_x],
                color=colors[i],
                marker=markers[j],
                label=label,
            )

            print(f"\n\n{label} got: err_x: {err_x}, err_u: {err_u}")

            del integrator

    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
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
    # plot_cstr(dt_plant, X_all, U_all, X_ref, U_ref, mpc_params.umin, mpc_params.umax, labels_all, x_min=[0, 250, 0.], x_max=[1, 500, 0.8])


if __name__ == "__main__":
    main()

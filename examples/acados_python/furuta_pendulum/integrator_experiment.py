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

import os, pickle
from furuta_model import get_furuta_model

import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from acados_template import latexify_plot, AcadosSim, AcadosSimSolver
from dataclasses import dataclass

latexify_plot()

X0 = np.array([0, np.pi/2, 0, 0])

WITH_FORWARD_SENSITIVITY = True
DT = 1e-3 * 60
# DT = 1e-3 * 2

@dataclass
class IntegratorSetting:
    integrator_type: str
    num_stages: int
    num_steps: int
    newton_iter: int
    newton_tol: float
    colloaction_type: str = "GAUSS_RADAU_IIA"
    jac_reuse: bool = True
    sens_forw: bool = True


def get_results_filename_from_setting(setting: IntegratorSetting):
    return os.path.join('results', f"dt_{DT}_{setting.integrator_type}_{setting.num_stages}_{setting.num_steps}_{setting.newton_iter}_{setting.newton_tol}{'_forw' if setting.sens_forw else ''}.pkl")


def get_order(setting: IntegratorSetting) -> int:
    if setting.integrator_type == "ERK":
        return setting.num_stages
    elif setting.integrator_type == "IRK":
        if setting.colloaction_type == "GAUSS_RADAU_IIA":
            return 2 * setting.num_stages - 1
        elif setting.colloaction_type == "GAUSS_LEGENDRE":
            return 2 * setting.num_stages

def setup_acados_integrator(model, dt, integrator_setting: IntegratorSetting):
    sim = AcadosSim()
    sim.model = model
    sim.solver_options.T = dt
    sim.solver_options.num_stages = integrator_setting.num_stages
    sim.solver_options.num_steps = integrator_setting.num_steps
    sim.solver_options.integrator_type = integrator_setting.integrator_type
    sim.solver_options.newton_iter = integrator_setting.newton_iter
    sim.solver_options.newton_tol = integrator_setting.newton_tol
    sim.solver_options.collocation_type = integrator_setting.colloaction_type
    sim.solver_options.sens_forw = True
    sim.solver_options.sens_adj = False
    sim.solver_options.sens_algebraic = False
    sim.solver_options.sens_hess = False
    sim.solver_options.sim_method_jac_reuse = integrator_setting.jac_reuse

    acados_integrator = AcadosSimSolver(sim, verbose=False)
    return acados_integrator


def simulate(
        integrator: AcadosSimSolver, x0, u_trajectory, n_runs=1
):
    Nsim = u_trajectory.shape[1]
    nx = integrator.acados_sim.dims.nx
    x_trajectory = np.zeros((nx, Nsim + 1))
    timings = 1e30 * np.ones((Nsim, n_runs))
    x_trajectory[:, 0] = x0

    for j in range(n_runs):
        for i in range(Nsim):
            x_trajectory[:, i + 1] = integrator.simulate(x_trajectory[:, i], u_trajectory[:, i], xdot=np.zeros((nx,)))
            timings[i, j] = integrator.get('time_tot')
    timings = np.mean(timings, axis=1)
    return x_trajectory, timings


def run_experiment(settings: list[IntegratorSetting]):
    dt_plant = DT

    model = get_furuta_model()
    Nsim = 1
    n_runs = 20

    ref_setting = IntegratorSetting(
        num_stages=8,
        num_steps=100,
        newton_iter=20,
        newton_tol=1e-16,
        integrator_type="IRK"
    )

    integrator = setup_acados_integrator(
        model,
        dt_plant,
        ref_setting
    )

    nu = 1
    np.random.seed(0)
    u_trajectory = 0.03 * (np.random.random((nu, Nsim)) - 0.5)
    X_exact, _ = simulate(integrator, X0, u_trajectory)
    del integrator

    for setting in settings:

        print(f"Running simulation with {setting}")

        model = get_furuta_model()

        integrator = setup_acados_integrator(
            model,
            dt_plant,
            setting
        )

        x_traj, timings_integrator = simulate(integrator, X0, u_trajectory, n_runs=n_runs)
        err_x = np.max(np.abs(x_traj - X_exact))

        # store results
        result = {"x_traj": x_traj, "timings_integrator": timings_integrator, "err_x": err_x}
        filename = get_results_filename_from_setting(setting)
        with open(filename, "wb") as f:
            pickle.dump(result, f)

        print(f"got: err_x: {err_x}")

        del integrator


def plot_results(settings: list[IntegratorSetting], results: list[dict]):

    SAVE_FIG = True
    plt.figure(figsize=(8, 5))

    integrator_type_vals = sorted(set([setting.integrator_type for setting in settings]))
    num_steps_vals = sorted(set([setting.num_steps for setting in settings]))
    order_vals = sorted(set([get_order(setting) for setting in settings]))

    markers = ["v", "s", "D", "^", ">", "<", "1", "2", "3", "4"]
    colors = [plt.cm.tab10(i) for i in range(len(order_vals))]
    fillstyles = ["full", "none"]

    for setting, result in zip(settings, results):
        timings_integrator = result["timings_integrator"]
        err_x = result["err_x"]
        order = get_order(setting)

        plt.plot(
            [np.mean(timings_integrator * 1e3)],
            [err_x],
            color=colors[order_vals.index(order)],
            marker=markers[num_steps_vals.index(setting.num_steps)],
            fillstyle=fillstyles[integrator_type_vals.index(setting.integrator_type)],
        )

    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
    plt.ylim([1e-16, 1e1])
    plt.ylabel(r"error $\Vert x - x^{\mathrm{exact}}\Vert_{\infty}$")
    plt.xlabel("mean computation time [ms]")

    legend_elements = \
        [
            Line2D([0], [0], color='k', lw=0, marker='o', label=f"{integrator_type}", fillstyle=fillstyles[i])
            for i, integrator_type in enumerate(integrator_type_vals)
        ] + \
        [
            Line2D(
                [0],
                [0],
                marker=markers[j],
                lw=0,
                color="k",
                label=f"{num_steps} steps",
            )
            for j, num_steps in enumerate(num_steps_vals)
        ] + \
        [
            Line2D([0], [0], color=colors[i], lw=4, label=f"order {order}")
            for i, order in enumerate(order_vals)
        ]

    plt.legend(handles=legend_elements, ncol=2)

    plt.title(f"Computation time vs. accuracy {'with' if WITH_FORWARD_SENSITIVITY else 'without'} sensitivity propagation, $\Delta t$ = {DT*1e3} ms")
    if SAVE_FIG:
        fig_filename = f"integrator_experiment_{DT}.pdf"
        plt.savefig(
            fig_filename, bbox_inches="tight", transparent=True, pad_inches=0.05
        )
        print(f"\nstored figure in {fig_filename}")

    plt.show()


def main():
    settings = []
    newton_iter = 20
    newton_tol = 1e-10
    for integrator_type in ["ERK", "IRK"]:
        if integrator_type == "ERK":
            num_stages_vals = [1, 2, 4]
            num_steps_vals = [1, 2, 5, 10, 20]
        elif integrator_type == "IRK":
            num_stages_vals = [1, 2, 3, 4]
            num_steps_vals = [1, 2, 5, 10, 20]
        for num_stages in num_stages_vals:
            for num_steps in num_steps_vals:
                settings.append(IntegratorSetting(num_stages=num_stages, num_steps=num_steps, integrator_type=integrator_type, newton_iter=newton_iter, newton_tol=newton_tol))

    run_experiment(settings)

    results = []
    for setting in settings:
        filename = get_results_filename_from_setting(setting)
        with open(filename, "rb") as f:
            results.append(pickle.load(f))
    plot_results(settings, results)

if __name__ == "__main__":
    main()

#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# authors: Anton Pozharskiy, Katrin Baumgaertner, Jonathan Frey

import os, pickle
from export_chain_mass_model import export_chain_mass_model
from utils import *

import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from acados_template import latexify_plot, AcadosSim, AcadosSimSolver
from dataclasses import dataclass

latexify_plot()

WITH_FORWARD_SENSITIVITY = True
N_MASSES = 10
#DT = 1e-1
DT = 2e-1
#DT = 1e-3 * 2

@dataclass
class IntegratorSetting:
    integrator_type: str
    num_stages: int
    num_steps: int
    newton_iter: int
    newton_tol: float
    collocation_type: str = "GAUSS_RADAU_IIA"
    jac_reuse: bool = False
    sens_forw: bool = WITH_FORWARD_SENSITIVITY


def get_results_filename_from_setting(setting: IntegratorSetting):
    if setting.integrator_type == "ERK":
        return os.path.join('results', f"dt_{DT}_{setting.integrator_type}_{setting.num_stages}_{setting.num_steps}_{setting.newton_iter}_{setting.newton_tol}{'_forw' if setting.sens_forw else ''}.pkl")
    else:
        return os.path.join('results', f"dt_{DT}_{setting.integrator_type}_{setting.collocation_type}_{setting.num_stages}_{setting.num_steps}_{setting.newton_iter}_{setting.newton_tol}{'_forw' if setting.sens_forw else ''}.pkl")


def get_order(setting: IntegratorSetting) -> int:
    if setting.integrator_type == "ERK":
        return setting.num_stages
    elif setting.integrator_type == "IRK":
        if setting.collocation_type == "GAUSS_RADAU_IIA":
            return 2 * setting.num_stages - 1
        elif setting.collocation_type == "GAUSS_LEGENDRE":
            return 2 * setting.num_stages
        elif setting.collocation_type == "SDIRK":
            return {1: 1, 2: 3, 3: 4, 9: 6}[setting.num_stages]
        elif setting.collocation_type == "DIRK":
            return {6: 6, 10: 7}[setting.num_stages]

def setup_acados_integrator(model, dt, integrator_setting: IntegratorSetting):
    sim = AcadosSim()
    sim.model = model
    sim.name = "test_example"
    sim.solver_options.T = dt
    sim.solver_options.num_stages = integrator_setting.num_stages
    sim.solver_options.num_steps = integrator_setting.num_steps
    sim.solver_options.integrator_type = integrator_setting.integrator_type
    sim.solver_options.newton_iter = integrator_setting.newton_iter
    sim.solver_options.newton_tol = integrator_setting.newton_tol
    sim.solver_options.collocation_type = integrator_setting.collocation_type
    sim.solver_options.sens_forw = integrator_setting.sens_forw
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
    S_forw = integrator.get("S_forw")
    return x_trajectory, S_forw, timings


def run_experiment(settings: list[IntegratorSetting]):
    dt_plant = DT

    chain_params = get_chain_params()
    chain_params["n_mass"] = N_MASSES
    n_mass = chain_params["n_mass"]
    M = chain_params["n_mass"] - 2 # number of intermediate masses
    # Ts = chain_params["Ts"]
    # Tsim = chain_params["Tsim"]
    N = chain_params["N"]
    # u_init = chain_params["u_init"]
    with_wall = chain_params["with_wall"]
    yPosWall = chain_params["yPosWall"]
    m = chain_params["m"]
    D = chain_params["D"]
    L = chain_params["L"]

    model = export_chain_mass_model(n_mass, m, D, L)
    Nsim = 1
    n_runs = 20

    nx = model.x.size1()
    nu = model.u.size1()

    # position of last mass
    xPosFirstMass = np.zeros((3,1))
    xEndRef = np.zeros((3,1))
    xEndRef[0] = L * (M+1) * 6
    # initial state
    pos0_x = np.linspace(xPosFirstMass[0], xEndRef[0], n_mass)
    x0 = np.zeros((nx,))
    x0[:3*(M+1):3] = pos0_x[1:].reshape((M+1))

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

    np.random.seed(0)
    u_trajectory = 0.03 * (np.random.random((nu, Nsim)) - 0.5)
    X_exact, S_forw_exact, _ = simulate(integrator, x0, u_trajectory)
    del integrator

    for setting in settings:

        print(f"Running simulation with {setting}")

        model = export_chain_mass_model(n_mass, m, D, L)

        integrator = setup_acados_integrator(
            model,
            dt_plant,
            setting
        )

        x_traj, S_forw, timings_integrator = simulate(integrator, x0, u_trajectory, n_runs=n_runs)
        err_x = np.max(np.abs(x_traj - X_exact))
        err_S_forw = np.max(np.abs(S_forw - S_forw_exact))

        # store results
        result = {"x_traj": x_traj, "timings_integrator": timings_integrator, "err_x": err_x, "err_S_forw": err_S_forw}
        filename = get_results_filename_from_setting(setting)
        with open(filename, "wb") as f:
            pickle.dump(result, f)

        print(f"got: err_x: {err_x}")
        print(f"got: err_S_forw: {err_S_forw}")

        del integrator


def plot_results(settings: list[IntegratorSetting], results: list[dict]):

    SAVE_FIG = True
    plt.figure(figsize=(8, 5))

    integrator_type_vals = sorted(set([(setting.integrator_type,setting.collocation_type) for setting in settings]))
    num_steps_vals = sorted(set([setting.num_steps for setting in settings]))
    order_vals = sorted(set([get_order(setting) for setting in settings]))

    markers = ["v", "s", "D", "^", ">", "<", "1", "2", "3", "4"]
    colors = [plt.cm.tab10(i) for i in range(len(order_vals))]
    fillstyles = ["full", "bottom", "none"]

    for setting, result in zip(settings, results):
        timings_integrator = result["timings_integrator"]
        err_x = result["err_x"]
        order = get_order(setting)

        plt.plot(
            [np.mean(timings_integrator * 1e3)],
            [err_x],
            color=colors[order_vals.index(order)],
            marker=markers[num_steps_vals.index(setting.num_steps)],
            fillstyle=fillstyles[integrator_type_vals.index((setting.integrator_type,setting.collocation_type))],
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
                label=f"{num_steps} steps" if num_steps > 1 else f"{num_steps} step",
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
    #for integrator_type in ["ERK", "IRK"]:
    for integrator_type in ["IRK"]:
        if integrator_type == "IRK":
            collocation_types = ["DIRK", "SDIRK", "GAUSS_RADAU_IIA"]
            #collocation_types = ["SDIRK", "GAUSS_RADAU_IIA"]
        else:
            collocation_types = ["GAUSS_RADAU_IIA"]
        for collocation_type in collocation_types:
            if integrator_type == "ERK":
                num_stages_vals = [1, 2, 4]
                num_steps_vals = [1, 2, 5, 10, 20]
            elif integrator_type == "IRK" and collocation_type == "SDIRK":
                num_stages_vals = [2, 3, 9]
                num_steps_vals = [1, 2, 5, 10, 20]
            elif integrator_type == "IRK" and collocation_type == "DIRK":
                num_stages_vals = [6, 10]
                num_steps_vals = [1, 2, 5, 10, 20]
            elif integrator_type == "IRK" and collocation_type == "GAUSS_RADAU_IIA":
                num_stages_vals = [2, 3, 4]
                num_steps_vals = [1, 2, 5, 10, 20]
            for num_stages in num_stages_vals:
                for num_steps in num_steps_vals:
                    settings.append(
                        IntegratorSetting(
                            num_stages=num_stages,
                            num_steps=num_steps,
                            integrator_type=integrator_type,
                            newton_iter=newton_iter,
                            newton_tol=newton_tol,
                            collocation_type=collocation_type,
                        )
                    )
    run_experiment(settings)

    results = []
    for setting in settings:
        filename = get_results_filename_from_setting(setting)
        with open(filename, "rb") as f:
            results.append(pickle.load(f))
    plot_results(settings, results)

if __name__ == "__main__":
    main()

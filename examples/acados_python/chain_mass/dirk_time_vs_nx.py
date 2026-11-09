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
from collections import defaultdict


latexify_plot()

WITH_FORWARD_SENSITIVITY = True
#DT = 1e-1
DT = 2e-1
#DT = 1e-3 * 2

@dataclass
class IntegratorSetting:
    N_masses: int
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
        return os.path.join('results', f"N_{setting.N_masses}_dt_{DT}_{setting.integrator_type}_{setting.num_stages}_{setting.num_steps}_{setting.newton_iter}_{setting.newton_tol}{'_forw' if setting.sens_forw else ''}.pkl")
    else:
        return os.path.join('results', f"N_{setting.N_masses}_dt_{DT}_{setting.integrator_type}_{setting.collocation_type}_{setting.num_stages}_{setting.num_steps}_{setting.newton_iter}_{setting.newton_tol}{'_forw' if setting.sens_forw else ''}.pkl")


def get_order(setting: IntegratorSetting) -> int:
    if setting.integrator_type == "ERK":
        return setting.num_stages
    elif setting.integrator_type == "IRK":
        if setting.collocation_type == "GAUSS_RADAU_IIA":
            return 2 * setting.num_stages - 1
        elif setting.collocation_type == "GAUSS_LEGENDRE":
            return 2 * setting.num_stages
        elif setting.collocation_type == "SDIRK":
            return {1: 1, 2: 3, 3: 4, 5: 5, 9: 6, 11: 7}[setting.num_stages]
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
    sim.solver_options.exploit_singly_implicit = True
    sim.solver_options.sim_method_jac_reuse = integrator_setting.jac_reuse

    acados_integrator = AcadosSimSolver(sim, verbose=False)
    return acados_integrator


def simulate(
        integrator: AcadosSimSolver, x0, u_trajectory, n_runs=1
):
    Nsim = u_trajectory.shape[1]
    nx = integrator.sim.dims.nx
    print(f"nx = {nx}")
    x_trajectory = np.zeros((nx, Nsim + 1))
    timings = 1e30 * np.ones((Nsim, n_runs))
    x_trajectory[:, 0] = x0

    for j in range(n_runs):
        print(f"run {j}")
        for i in range(Nsim):
            x_trajectory[:, i + 1] = integrator.simulate(x_trajectory[:, i], u_trajectory[:, i], xdot=np.zeros((nx,)))
            timings[i, j] = integrator.get('time_tot')
    timings = np.median(timings, axis=1)
    S_forw = integrator.get("S_forw")
    return x_trajectory, S_forw, timings


def run_experiment(settings: list[IntegratorSetting]):
    dt_plant = DT
    for setting in settings:
        chain_params = get_chain_params()
        chain_params["n_mass"] = setting.N_masses
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
        n_runs = 10

        nx = model.x.size1()
        nu = model.u.size1()

        np.random.seed(0)
        u_trajectory = 0.03 * (np.random.random((nu, Nsim)) - 0.5)

        # position of last mass
        xPosFirstMass = np.zeros((3,1))
        xEndRef = np.zeros((3,1))
        xEndRef[0] = L * (M+1) * 6
        # initial state
        pos0_x = np.linspace(xPosFirstMass[0], xEndRef[0], n_mass)
        x0 = np.zeros((nx,))
        x0[:3*(M+1):3] = pos0_x[1:].reshape((M+1))

        print(f"Running simulation with {setting}")

        model = export_chain_mass_model(n_mass, m, D, L)

        integrator = setup_acados_integrator(
            model,
            dt_plant,
            setting
        )

        x_traj, S_forw, timings_integrator = simulate(integrator, x0, u_trajectory, n_runs=n_runs)

        # store results
        result = {"timings_integrator": timings_integrator, "nx": nx}
        filename = get_results_filename_from_setting(setting)
        with open(filename, "wb") as f:
            pickle.dump(result, f)

        del integrator


def plot_results(settings: list[IntegratorSetting], results: list[dict]):

    SAVE_FIG = True
    plt.figure(figsize=(5, 4))

    integrator_type_vals = sorted(set([(setting.integrator_type,setting.collocation_type) for setting in settings]))
    num_steps_vals = sorted(set([setting.num_steps for setting in settings]))
    order_vals = sorted(set([get_order(setting) for setting in settings]))
    nx_vals = sorted(set([result["nx"] for result in results]))

    markers = ["v", "s", "D", "^", ">", "<", "1", "2", "3", "4"]
    colors = [plt.cm.tab10(i) for i in range(len(order_vals))]
    fillstyles = ["full","none","bottom"]


    res_vals = defaultdict(list)
    for setting,result in zip(settings,results):
        timings_integrator = result["timings_integrator"]
        res_vals[(setting.num_stages, setting.collocation_type)].append(timings_integrator)

    ratio_vals = defaultdict(np.array)

    # get ratio for
    ratio_vals["Order 3"] = np.sort(res_vals[(2,"SDIRK")]) / np.sort(res_vals[(2,"GAUSS_RADAU_IIA")])
    ratio_vals["Order 4"] = np.sort(res_vals[(3,"SDIRK")]) / np.sort(res_vals[(2,"GAUSS_LEGENDRE")])
    ratio_vals["Order 5"] = np.sort(res_vals[(5,"SDIRK")]) / np.sort(res_vals[(3,"GAUSS_RADAU_IIA")])
    ratio_vals["Order 6"] = np.sort(res_vals[(9,"SDIRK")]) / np.sort(res_vals[(3,"GAUSS_LEGENDRE")])
    ratio_vals["Order 7"] = np.sort(res_vals[(11,"SDIRK")]) / np.sort(res_vals[(4,"GAUSS_RADAU_IIA")])

    ic = 0
    for k,v in ratio_vals.items():
        plt.plot(nx_vals, v, label=str(k), color=f"C{ic}", marker="o")
        ic += 1

    plt.grid()
    #plt.yscale("log")
    plt.xscale("log")
    plt.ylim([0, 1.5])
    plt.ylabel(r"computation time ratio")
    plt.xlabel("$n_x$")

    # legend_elements = \
    #     [
    #         Line2D([0], [0], color='k', lw=0, marker='o', label=f"{integrator_type}", fillstyle=fillstyles[i])
    #         for i, integrator_type in enumerate(integrator_type_vals)
    #     ] + \
    #     [
    #         Line2D(
    #             [0],
    #             [0],
    #             marker=markers[j],
    #             lw=0,
    #             color="k",
    #             label=f"{num_steps} steps" if num_steps > 1 else f"{num_steps} step",
    #         )
    #         for j, num_steps in enumerate(num_steps_vals)
    #     ] + \
    #     [
    #         Line2D([0], [0], color=colors[i], lw=4, label=f"order {order}")
    #         for i, order in enumerate(order_vals)
    #     ]

    plt.legend()

    #plt.title(f"Computation time vs. State Size")
    ax = plt.gca()
    ax.set_xticklabels([], minor=True)
    if SAVE_FIG:
        fig_filename = f"t_vs_dx.pdf"
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
    for N_masses in [5, 7, 10, 25, 50, 75, 100, 150]:
        for (n_stages,collocation_type) in [(2,"SDIRK"), (3,"SDIRK"), (5,"SDIRK"), (9,"SDIRK"), (11, "SDIRK"), (2,"GAUSS_RADAU_IIA"), (3,"GAUSS_RADAU_IIA"),  (4,"GAUSS_RADAU_IIA"), (2,"GAUSS_LEGENDRE"), (3,"GAUSS_LEGENDRE"),]:
             settings.append(
                 IntegratorSetting(
                     N_masses=N_masses,
                     num_stages=n_stages,
                     num_steps=3,
                     integrator_type="IRK",
                     newton_iter=newton_iter,
                     newton_tol=newton_tol,
                     collocation_type=collocation_type,
                 )
             )
    #run_experiment(settings)

    results = []
    for setting in settings:
        filename = get_results_filename_from_setting(setting)
        with open(filename, "rb") as f:
            results.append(pickle.load(f))
    plot_results(settings, results)

if __name__ == "__main__":
    main()

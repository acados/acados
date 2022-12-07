import sys
import os
from xmlrpc.client import Boolean
local_path = os.path.dirname(os.path.abspath(__file__))
cstr_source_dir = os.path.join(local_path, '..', '..', 'cstr')
sys.path.append(cstr_source_dir)
zoro_source_dir = os.path.join(local_path, '..')
sys.path.append(zoro_source_dir)

import numpy as np
import matplotlib
from typing import Optional

from cstr_model import CSTRParameters, setup_cstr_model, setup_linearized_model
from zoro_ocp_solver import MpcCSTRParameters, setup_acados_ocp_solver, AcadosOcpSolver
from setup_acados_integrator import setup_acados_integrator, AcadosSimSolver
from cstr_utils import plot_cstr, compute_lqr_gain


def simulate(
    controller: Optional[AcadosOcpSolver],
    plant: AcadosSimSolver,
    x0: np.ndarray,
    Nsim: int,
    X_ref: np.ndarray,
    U_ref: np.ndarray,
    cstr_tightening: Boolean
):

    nx = X_ref.shape[1]
    nu = U_ref.shape[1]

    X = np.ndarray((Nsim + 1, nx))
    U = np.ndarray((Nsim, nu))
    timings_solver = np.zeros((Nsim))
    timings_integrator = np.zeros((Nsim))

    # closed loop
    xcurrent = x0
    X[0, :] = xcurrent

    for i in range(Nsim):

        # set initial state
        controller.set(0, "lbx", xcurrent)
        controller.set(0, "ubx", xcurrent)

        yref = np.concatenate((X_ref[i, :], U_ref[i, :]))
        for stage in range(controller.acados_ocp.dims.N):
            controller.set(stage, "yref", yref)
        controller.set(controller.acados_ocp.dims.N, "yref", X_ref[i, :])

        if cstr_tightening:
            controller.custom_update([])

        # solve ocp
        status = controller.solve()

        if status != 0:
            controller.print_statistics()
            raise Exception(
                f"acados controller returned status {status} in simulation step {i}. Exiting."
            )

        U[i, :] = controller.get(0, "u")
        timings_solver[i] = controller.get_stats("time_tot")

        # simulate system
        plant.set("x", xcurrent)
        plant.set("u", U[i, :])

        if plant.acados_sim.solver_options.integrator_type == "IRK":
            plant.set("xdot", np.zeros((nx,)))

        status = plant.solve()
        if status != 0:
            raise Exception(
                f"acados integrator returned status {status} in simulation step {i}. Exiting."
            )

        timings_integrator[i] = plant.get("time_tot")
        # update state
        xcurrent = plant.get("x")
        X[i + 1, :] = xcurrent

    return X, U, timings_solver, timings_integrator



def main():

    Tsim = 25
    dt_plant = 0.25  # [min]

    Nsim = int(Tsim / dt_plant)
    if not (Tsim / dt_plant).is_integer():
        print("WARNING: Tsim / dt_plant should be an integer!")

    # steady-state
    xs = np.array([[0.878, 324.5, 0.659]]).T
    us = np.array([[300, 0.1]]).T
    # initial state
    x0 = np.array([0.05, 0.75, 0.5]) * xs.ravel()
    # constant ref
    X_ref = np.tile(xs, Nsim + 1).T
    U_ref = np.tile(us, Nsim).T

    cstr_params = CSTRParameters()
    mpc_params = MpcCSTRParameters(xs=cstr_params.xs, us=cstr_params.us)
    model = setup_cstr_model(cstr_params)
    plant_model = setup_cstr_model(cstr_params)

    integrator = setup_acados_integrator(plant_model, dt_plant, cstr_param=cstr_params)

    X_all = []
    U_all = []
    labels_all = []
    timings_solver_all = []

    # simulation with Nominal MPC
    label = "Nominal"
    print(f"\n\nRunning simulation with {label}\n\n")
    mpc_params.cstr_tightening = False
    ocp_solver = setup_acados_ocp_solver(model, mpc_params, cstr_params=cstr_params)

    X, U, timings_solver, _ = simulate(
        ocp_solver, integrator, x0, Nsim, X_ref=X_ref, U_ref=U_ref, cstr_tightening=False
    )
    X_all.append(X)
    U_all.append(U)
    timings_solver_all.append(timings_solver)
    labels_all.append(label)
    ocp_solver = None

    # simulation with zoRO MPC
    label = "zoRO"
    print(f"\n\nRunning simulation with {label}\n\n")
    mpc_params.cstr_tightening = True
    ocp_solver = setup_acados_ocp_solver(model, mpc_params, cstr_params=cstr_params)

    X, U, timings_solver, _ = simulate(
        ocp_solver, integrator, x0, Nsim, X_ref=X_ref, U_ref=U_ref, cstr_tightening=True
    )
    X_all.append(X)
    U_all.append(U)
    timings_solver_all.append(timings_solver)
    labels_all.append(label)
    ocp_solver = None

    # Evaluation
    print("\nTiming evaluation:\n------------------")
    for i in range(len(labels_all)):
        label = labels_all[i]
        timings_solver = timings_solver_all[i] * 1e3
        print(
            f"{label}:\n min: {np.min(timings_solver):.3f} ms, mean: {np.mean(timings_solver):.3f} ms, max: {np.max(timings_solver):.3f} ms\n"
        )

    matplotlib.use("TkAgg")
    plot_cstr(
        dt_plant,
        X_all,
        U_all,
        X_ref,
        U_ref,
        mpc_params.umin,
        mpc_params.umax,
        labels_all,
    )  # , fig_filename='cstr_acados_RTI.pdf')


if __name__ == "__main__":
    main()
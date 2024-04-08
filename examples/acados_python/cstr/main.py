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

from cstr_model import CstrParameters, setup_cstr_model, setup_linearized_model
from setup_acados_ocp_solver import (
    MpcCstrParameters,
    setup_acados_ocp_solver,
    AcadosOcpSolver,
)
from setup_acados_integrator import setup_acados_integrator, AcadosSimSolver
import numpy as np
import casadi as ca
from cstr_utils import plot_cstr
from typing import Optional


def simulate(
    controller: Optional[AcadosOcpSolver],
    plant: AcadosSimSolver,
    x0: np.ndarray,
    Nsim: int,
    X_ref: np.ndarray,
    U_ref: np.ndarray,
    with_reference_profile: bool = False,
):

    nx = X_ref.shape[1]
    nu = U_ref.shape[1]

    X = np.zeros((Nsim + 1, nx))
    U = np.zeros((Nsim, nu))
    timings_solver = np.zeros((Nsim))
    timings_integrator = np.zeros((Nsim))

    # closed loop
    xcurrent = x0
    X[0, :] = xcurrent

    for i in range(Nsim):

        if controller is None:
            U[i, :] = U_ref[i, :]
        else:
            if with_reference_profile:
                t0 = i * plant.acados_sim.solver_options.T
                for stage, dt in enumerate(controller.acados_ocp.solver_options.shooting_nodes):
                    t = t0 + dt
                    controller.set_params_sparse(stage, np.array([1]), np.array([t]))
            else:
                yref = np.concatenate((X_ref[i, :], U_ref[i, :]))
                for stage in range(controller.acados_ocp.dims.N):
                    controller.set(stage, "yref", yref)
                controller.set(controller.acados_ocp.dims.N, "yref", X_ref[i, :])

            # solve ocp
            U[i, :] = controller.solve_for_x0(xcurrent)
            timings_solver[i] = controller.get_stats("time_tot")

        # simulate system
        xcurrent = plant.simulate(x=xcurrent, u=U[i, :], xdot=np.zeros((nx,)))
        timings_integrator[i] = plant.get("time_tot")
        X[i + 1, :] = xcurrent

    return X, U, timings_solver, timings_integrator


def main():
    with_nmpc = True
    with_timevar_ref_nmpc = True
    with_linear_mpc = True
    with_nmpc_rti = True

    Tsim = 25
    dt_plant = 0.25  # [min]

    cstr_params = CstrParameters()
    mpc_params = MpcCstrParameters(xs=cstr_params.xs, us=cstr_params.us)
    model = setup_cstr_model(cstr_params)
    linearized_model = setup_linearized_model(model, cstr_params, mpc_params)
    plant_model = setup_cstr_model(cstr_params)

    Nsim = int(Tsim / dt_plant)
    if not (Tsim / dt_plant).is_integer():
        print("WARNING: Tsim / dt_plant should be an integer!")

    integrator = setup_acados_integrator(plant_model, dt_plant, cstr_param=cstr_params)

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

    X_all = []
    U_all = []
    labels_all = []
    timings_solver_all = []

    # simulation with constant reference input
    label = "constant reference input"
    print(f"\n\nRunning simulation with {label}\n\n")
    X, U, timings_solver, _ = simulate(None, integrator, x0, Nsim, X_ref, U_ref)
    X_all.append(X)
    U_all.append(U)
    timings_solver_all.append(timings_solver)
    labels_all.append(label)

    # simulation with NMPC controller
    if with_nmpc:
        label = "NMPC"
        print(f"\n\nRunning simulation with {label}\n\n")
        ocp_solver = setup_acados_ocp_solver(model, mpc_params, cstr_params=cstr_params)

        X, U, timings_solver, _ = simulate(
            ocp_solver, integrator, x0, Nsim, X_ref=X_ref, U_ref=U_ref
        )
        X_all.append(X)
        U_all.append(U)
        timings_solver_all.append(timings_solver)
        labels_all.append(label)
        ocp_solver = None

    # simulation with time varying reference NMPC controller
    if with_timevar_ref_nmpc:
        label = "NMPC with time-varying reference"
        ocp_model = setup_cstr_model(cstr_params)
        ocp_model.t = ca.SX.sym("t")
        t0 = ca.SX.sym("t0")
        ocp_model.p = ca.vertcat(ocp_model.p, t0)
        t = t0 + ocp_model.t

        tjump1 = Njump * dt_plant
        tjump2 = 2 * Njump * dt_plant
        reference_profile = ca.if_else(t < tjump1, ca.vertcat(xs, us), ca.vertcat(xs2, us2))
        reference_profile = ca.if_else(t < tjump2, reference_profile, ca.vertcat(xs, us))
        print(f"\n\nRunning simulation with {label}\n\n")

        ocp_solver = setup_acados_ocp_solver(ocp_model, mpc_params, cstr_params=cstr_params, reference_profile=reference_profile, cost_integration=True)

        X, U, timings_solver, _ = simulate(
            ocp_solver, integrator, x0, Nsim, X_ref=X_ref, U_ref=U_ref, with_reference_profile=True
        )
        X_all.append(X)
        U_all.append(U)
        timings_solver_all.append(timings_solver)
        labels_all.append(label)
        ocp_solver = None

    # simulation with LMPC controller
    if with_linear_mpc:
        label = "LMPC"
        print(f"\n\nRunning simulation with {label}\n\n")
        mpc_params.linear_mpc = True
        ocp_solver = setup_acados_ocp_solver(
            linearized_model, mpc_params, cstr_params=cstr_params, use_rti=True
        )
        mpc_params.linear_mpc = False

        X, U, timings_solver, _ = simulate(
            ocp_solver, integrator, x0, Nsim, X_ref=X_ref, U_ref=U_ref
        )
        X_all.append(X)
        U_all.append(U)
        timings_solver_all.append(timings_solver)
        labels_all.append(label)
        ocp_solver = None

    # simulation with NMPC RTI controller
    if with_nmpc_rti:
        label = "NMPC-RTI"
        print(f"\n\nRunning simulation with {label}\n\n")
        ocp_solver = setup_acados_ocp_solver(
            model, mpc_params, cstr_params=cstr_params, use_rti=True
        )

        X, U, timings_solver, _ = simulate(
            ocp_solver, integrator, x0, Nsim, X_ref=X_ref, U_ref=U_ref
        )
        X_all.append(X)
        U_all.append(U)
        timings_solver_all.append(timings_solver)
        labels_all.append(label)
        ocp_solver = None

    # Evaluation
    max_label_length = max([len(l) for l in labels_all])
    print(f"\n{'Timings in ms':{max_label_length}}  | {'miniumum':<8} | {'mean':<8} | {'maximum':<8} \n------------------")

    for i in range(len(labels_all)):
        label = labels_all[i]
        timings_solver = timings_solver_all[i] * 1e3  # scale to milliseconds
        print(
            f"{label:{max_label_length}}   {f'{np.min(timings_solver):.3f}':>10} {f'{np.mean(timings_solver):.3f}':>10} {f'{np.max(timings_solver):.3f}':>10}"
        )

    # plot results
    plot_cstr(
        dt_plant,
        X_all,
        U_all,
        X_ref,
        U_ref,
        mpc_params.umin,
        mpc_params.umax,
        labels_all,
        fig_filename='cstr_acados_RTI.pdf',
    )


if __name__ == "__main__":
    main()

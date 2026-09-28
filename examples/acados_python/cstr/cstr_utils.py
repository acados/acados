#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import matplotlib.pyplot as plt
import numpy as np
from acados_template import latexify_plot

latexify_plot()

def plot_cstr(
    dt,
    X_list,
    U_list,
    X_ref,
    U_ref,
    u_min,
    u_max,
    labels_list,
    fig_filename=None,
    x_min=None,
    x_max=None,
):
    nx = X_list[0].shape[1]
    nu = U_list[0].shape[1]

    Nsim = U_list[0].shape[0]

    ts = dt * np.arange(0, Nsim + 1)

    x_labels = ["$c$ [kmol/m$^3$]", "$T$ [K]", "$h$ [m]"]
    controls_lables = ["$T_c$ [K]", "$F$ [m$^3$/min]"]

    fig, axes = plt.subplots(ncols=2, nrows=nx, figsize=(10, 7))

    for i in range(nx):
        for X, label in zip(X_list, labels_list):
            axes[i, 0].plot(ts, X[:, i], label=label, alpha=0.7)

        axes[i, 0].step(
            ts,
            X_ref[:, i],
            alpha=0.8,
            where="post",
            label="reference",
            linestyle="dotted",
            color="k",
        )
        axes[i, 0].set_ylabel(x_labels[i])
        axes[i, 0].grid()
        axes[i, 0].set_xlim(ts[0], ts[-1])

        if x_max is not None:
            axes[i, 0].hlines(
                x_max[i], ts[0], ts[-1], linestyles="dashed", alpha=0.8, color="k"
            )

        if x_min is not None:
            axes[i, 0].set_ylim(bottom=x_min[i])

        if x_max is not None:
            axes[i, 0].set_ylim(top=x_max[i] * 1.05)

    for i in range(nu):
        for U, label in zip(U_list, labels_list):
            axes[i, 1].step(ts, np.append([U[0, i]], U[:, i]), label=label, alpha=0.7)
        axes[i, 1].step(
            ts,
            np.append([U_ref[0, i]], U_ref[:, i]),
            alpha=0.8,
            label="reference",
            linestyle="dotted",
            color="k",
        )
        axes[i, 1].set_ylabel(controls_lables[i])
        axes[i, 1].grid()

        axes[i, 1].hlines(
            u_max[i], ts[0], ts[-1], linestyles="dashed", alpha=0.8, color="k"
        )
        axes[i, 1].hlines(
            u_min[i], ts[0], ts[-1], linestyles="dashed", alpha=0.8, color="k"
        )
        axes[i, 1].set_xlim(ts[0], ts[-1])
        axes[i, 1].set_ylim(bottom=0.98 * u_min[i], top=1.02 * u_max[i])

    axes[1, 1].legend(bbox_to_anchor=(0.5, -1.7), loc="lower center")
    axes[-1, 0].set_xlabel("$t$ [min]")
    axes[1, 1].set_xlabel("$t$ [min]")

    fig.delaxes(axes[-1, 1])

    plt.subplots_adjust(
        left=None, bottom=None, right=None, top=None, hspace=0.3, wspace=0.4
    )
    if fig_filename is not None:
        plt.savefig(
            fig_filename, bbox_inches="tight", transparent=True, pad_inches=0.05
        )
        print(f"\nstored figure in {fig_filename}")

    plt.show()


def compute_lqr_gain(model, cstr_params, mpc_params):
    from scipy.linalg import solve_discrete_are

    # linearize dynamics
    from setup_acados_integrator import setup_acados_integrator

    integrator = setup_acados_integrator(
        model, mpc_params.dt, cstr_param=cstr_params, sensitivity_propagation=True
    )
    integrator.set("x", cstr_params.xs)
    integrator.set("u", cstr_params.us)

    integrator.solve()

    A_mat = integrator.get("Sx")
    B_mat = integrator.get("Su")

    Q_mat = mpc_params.dt * mpc_params.Q
    R_mat = mpc_params.dt * mpc_params.R
    P_mat = solve_discrete_are(A_mat, B_mat, Q_mat, R_mat)

    return A_mat, B_mat, P_mat

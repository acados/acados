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


"""
Test for solution sensitivities with many parameters.
"""

import os
import numpy as np
import casadi as ca
from casadi import SX, norm_2, vertcat
from casadi.tools import struct_symSX, entry
from casadi.tools.structure3 import DMStruct
import matplotlib.pyplot as plt
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
from utils import get_chain_params
from typing import Tuple
from plot_utils import plot_timings


def export_discrete_erk4_integrator_step(f_expl: SX, x: SX, u: SX, p: struct_symSX, h: float, n_stages: int = 2) -> ca.SX:
    """Define ERK4 integrator for continuous dynamics."""
    dt = h / n_stages
    ode = ca.Function("f", [x, u, p], [f_expl])
    xnext = x
    for _ in range(n_stages):
        k1 = ode(xnext, u, p)
        k2 = ode(xnext + dt / 2 * k1, u, p)
        k3 = ode(xnext + dt / 2 * k2, u, p)
        k4 = ode(xnext + dt * k3, u, p)
        xnext = xnext + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    return xnext


def define_param_struct_symSX(n_mass: int, disturbance: bool = True) -> DMStruct:
    """Define parameter struct."""
    n_link = n_mass - 1

    nx, nu = define_nx_nu(n_mass)

    param_entries = [
        entry("m", shape=1, repeat=n_link),
        entry("D", shape=3, repeat=n_link),
        entry("L", shape=3, repeat=n_link),
        entry("C", shape=3, repeat=n_link),
        entry("Q", shape=(nx, nx)),
        entry("R", shape=(nu, nu)),
    ]

    if disturbance:
        param_entries.append(entry("w", shape=3, repeat=n_mass - 2))

    return struct_symSX(param_entries)


def define_nx_nu(n_mass: int) -> Tuple[int, int]:
    """Define number of states and control inputs."""
    M = n_mass - 2  # number of intermediate masses
    nx = (2 * M + 1) * 3  # differential states
    nu = 3  # control inputs

    return nx, nu


def find_idx_for_labels(sub_vars: SX, sub_label: str) -> list[int]:
    """Return a list of indices where sub_label is part of the variable label."""
    return [i for i, label in enumerate(sub_vars.str().strip("[]").split(", ")) if sub_label in label]


def export_chain_mass_model(n_mass: int, Ts: float = 0.2, disturbance: bool = False) -> Tuple[AcadosModel, DMStruct]:
    """Export chain mass model for acados."""
    x0 = np.array([0, 0, 0])  # fix mass (at wall)

    M = n_mass - 2  # number of intermediate massesu

    nx, nu = define_nx_nu(n_mass)

    xpos = SX.sym("xpos", (M + 1) * 3, 1)  # position of fixed mass eliminated
    xvel = SX.sym("xvel", M * 3, 1)
    u = SX.sym("u", nu, 1)
    xdot = SX.sym("xdot", nx, 1)

    f = SX.zeros(3 * M, 1)  # force on intermediate masses
    p = define_param_struct_symSX(n_mass=n_mass, disturbance=disturbance)

    # Gravity force
    for i in range(M):
        f[3 * i + 2] = -9.81

    # Spring force
    for i in range(M + 1):
        if i == 0:
            dist = xpos[i * 3 : (i + 1) * 3] - x0
        else:
            dist = xpos[i * 3 : (i + 1) * 3] - xpos[(i - 1) * 3 : i * 3]

        F = ca.SX.zeros(3, 1)
        for j in range(3):
            F[j] = p["D", i, j] / p["m", i] * (1 - p["L", i, j] / norm_2(dist)) * dist[j]

        # mass on the right
        if i < M:
            f[i * 3 : (i + 1) * 3] -= F

        # mass on the left
        if i > 0:
            f[(i - 1) * 3 : i * 3] += F

    # Damping force
    for i in range(M + 1):
        if i == 0:
            vel = xvel[i * 3 : (i + 1) * 3]
        elif i == M:
            vel = u - xvel[(i - 1) * 3 : i * 3]
        else:
            vel = xvel[i * 3 : (i + 1) * 3] - xvel[(i - 1) * 3 : i * 3]

        F = ca.SX.zeros(3, 1)
        for j in range(3):
            F[j] = p["C", i, j] * vel[j]

        # mass on the right
        if i < M:
            f[i * 3 : (i + 1) * 3] -= F

        # mass on the left
        if i > 0:
            f[(i - 1) * 3 : i * 3] += F

    x = vertcat(xpos, xvel)

    # dynamics
    if disturbance:
        # Disturbance on intermediate masses
        for i in range(M):
            f[i * 3 : (i + 1) * 3] += p["w", i]
        model_name = "chain_mass_ds_" + str(n_mass)
    else:
        model_name = "chain_mass_" + str(n_mass)

    f_expl = vertcat(xvel, u, f)
    f_impl = xdot - f_expl
    f_disc = export_discrete_erk4_integrator_step(f_expl, x, u, p, Ts)

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.disc_dyn_expr = f_disc
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p.cat
    model.name = model_name

    p_map = p(0)

    return model, p_map


def compute_parametric_steady_state(
    model: AcadosModel, p: DMStruct, xPosFirstMass: np.ndarray, xEndRef: np.ndarray
) -> np.ndarray:
    """Compute steady state for chain mass model."""
    # TODO reuse/adapt the compute_steady_state function in utils.py

    p_ = p(0)
    p_["m"] = p["m"]
    p_["D"] = p["D"]
    p_["L"] = p["L"]
    p_["C"] = p["C"]

    nx = model.x.shape[0]

    # Free masses
    M = int((nx / 3 - 1) / 2)

    # initial guess for state
    pos0_x = np.linspace(xPosFirstMass[0], xEndRef[0], M + 2)
    x0 = np.zeros((nx, 1))
    x0[: 3 * (M + 1) : 3] = pos0_x[1:].reshape((M + 1, 1))

    # decision variables
    w = [model.x, model.xdot, model.u]

    # initial guess
    w0 = ca.vertcat(*[x0, np.zeros(model.xdot.shape), np.zeros(model.u.shape)])

    # constraints
    g = []
    g += [model.f_impl_expr]  # steady state
    g += [model.x[3 * M : 3 * (M + 1)] - xEndRef]  # fix position of last mass
    g += [model.u]  # don't actuate controlled mass

    # misuse IPOPT as nonlinear equation solver
    nlp = {"x": ca.vertcat(*w), "f": 0, "g": ca.vertcat(*g), "p": model.p}

    solver = ca.nlpsol("solver", "ipopt", nlp)
    sol = solver(x0=w0, lbg=0, ubg=0, p=p_.cat)

    x_ss = sol["x"].full()[:nx]

    return x_ss


def export_parametric_ocp(
    chain_params_: dict,
    qp_solver_ric_alg: int = 0,
    qp_solver: str = "PARTIAL_CONDENSING_HPIPM",
    hessian_approx: str = "GAUSS_NEWTON",
    integrator_type: str = "IRK",
    nlp_solver_type: str = "SQP",
    nlp_iter: int = 50,
    nlp_tol: float = 1e-5,
    random_scale: dict = {"m": 0.0, "D": 0.0, "L": 0.0, "C": 0.0},
) -> Tuple[AcadosOcp, DMStruct]:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()
    ocp.dims.N = chain_params_["N"]

    # export model
    ocp.model, p = export_chain_mass_model(n_mass=chain_params_["n_mass"], Ts=chain_params_["Ts"], disturbance=True)

    # parameters
    np.random.seed(chain_params_["seed"])
    m = [
        np.random.normal(chain_params_["m"], random_scale["m"] * chain_params_["m"])
        for _ in range(chain_params_["n_mass"] - 1)
    ]
    D = [
        np.random.normal(chain_params_["D"], random_scale["D"] * chain_params_["D"], 3)
        for _ in range(chain_params_["n_mass"] - 1)
    ]
    L = [
        np.random.normal(chain_params_["L"], random_scale["L"] * chain_params_["L"], 3)
        for _ in range(chain_params_["n_mass"] - 1)
    ]
    C = [
        np.random.normal(chain_params_["C"], random_scale["C"] * chain_params_["C"], 3)
        for _ in range(chain_params_["n_mass"] - 1)
    ]

    # Inermediate masses
    for i_mass in range(chain_params_["n_mass"] - 1):
        p["m", i_mass] = m[i_mass]
    for i_mass in range(chain_params_["n_mass"] - 2):
        p["w", i_mass] = np.zeros(3)

    # Links
    for i_link in range(chain_params_["n_mass"] - 1):
        p["D", i_link] = D[i_link]
        p["L", i_link] = L[i_link]
        p["C", i_link] = C[i_link]

    # initial state
    x_ss = compute_parametric_steady_state(
        model=ocp.model,
        p=p,
        xPosFirstMass=chain_params_["xPosFirstMass"],
        xEndRef=np.array([chain_params_["L"] * (chain_params_["n_mass"] - 1) * 6, 0.0, 0.0]),
    )
    x0 = x_ss

    nx = ocp.model.x.rows()
    nu = ocp.model.u.rows()

    M = chain_params_["n_mass"] - 2  # number of intermediate masses
    ocp.cost.cost_type = "EXTERNAL"
    ocp.cost.cost_type_e = "EXTERNAL"

    x_e = ocp.model.x - x_ss
    u_e = ocp.model.u - np.zeros((nu, 1))

    idx = find_idx_for_labels(define_param_struct_symSX(chain_params_["n_mass"], disturbance=True).cat, "Q")
    Q_sym = ca.reshape(ocp.model.p[idx], (nx, nx))
    q_diag = np.ones((nx, 1))
    q_diag[3 * M : 3 * M + 3] = M + 1
    p["Q"] = 2 * np.diagflat(q_diag)

    idx = find_idx_for_labels(define_param_struct_symSX(chain_params_["n_mass"], disturbance=True).cat, "R")
    R_sym = ca.reshape(ocp.model.p[idx], (nu, nu))
    p["R"] = 2 * np.diagflat(1e-2 * np.ones((nu, 1)))

    ocp.model.cost_expr_ext_cost = 0.5 * (x_e.T @ Q_sym @ x_e + u_e.T @ R_sym @ u_e)
    ocp.model.cost_expr_ext_cost_e = 0.5 * (x_e.T @ Q_sym @ x_e)

    ocp.model.cost_y_expr = vertcat(x_e, u_e)

    ocp.parameter_values = p.cat.full().flatten()

    # set constraints
    umax = 1 * np.ones((nu,))

    ocp.constraints.lbu = -umax
    ocp.constraints.ubu = umax
    ocp.constraints.x0 = x0.reshape((nx,))
    ocp.constraints.idxbu = np.array(range(nu))

    # solver options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = hessian_approx
    ocp.solver_options.integrator_type = integrator_type
    ocp.solver_options.nlp_solver_type = nlp_solver_type
    ocp.solver_options.sim_method_num_stages = 2
    ocp.solver_options.sim_method_num_steps = 2
    ocp.solver_options.nlp_solver_max_iter = nlp_iter

    if hessian_approx == "EXACT":
        ocp.solver_options.nlp_solver_step_length = 0.0
        ocp.solver_options.nlp_solver_max_iter = 1
        ocp.solver_options.qp_solver_iter_max = 200
        ocp.solver_options.tol = 1e-10
        ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg
        ocp.solver_options.qp_solver_cond_N = ocp.dims.N
        ocp.solver_options.with_solution_sens_wrt_params = True
    else:
        ocp.solver_options.nlp_solver_max_iter = nlp_iter
        ocp.solver_options.qp_solver_cond_N = ocp.dims.N
        ocp.solver_options.qp_tol = nlp_tol
        ocp.solver_options.tol = nlp_tol

    ocp.solver_options.tf = ocp.dims.N * chain_params_["Ts"]

    return ocp, p


def main_parametric(qp_solver_ric_alg: int = 0, chain_params_: dict = get_chain_params(), generate_code: bool = True) -> None:
    ocp, parameter_values = export_parametric_ocp(
        chain_params_=chain_params_, qp_solver_ric_alg=qp_solver_ric_alg, integrator_type="DISCRETE",
    )

    ocp_json_file = "acados_ocp_" + ocp.model.name + ".json"

    # Check if json_file exists
    if not generate_code and os.path.exists(ocp_json_file):
        ocp_solver = AcadosOcpSolver(ocp, json_file=ocp_json_file, build=False, generate=False)
    else:
        ocp_solver = AcadosOcpSolver(ocp, json_file=ocp_json_file)

    sensitivity_ocp, _ = export_parametric_ocp(
        chain_params_=chain_params_,
        qp_solver_ric_alg=qp_solver_ric_alg,
        hessian_approx="EXACT",
        integrator_type="DISCRETE",
    )
    sensitivity_ocp.model.name = f"{ocp.model.name}_sensitivity"

    ocp_json_file = "acados_sensitivity_ocp_" + sensitivity_ocp.model.name + ".json"
    # Check if json_file exists
    if not generate_code and os.path.exists(ocp_json_file):
        sensitivity_solver = AcadosOcpSolver(sensitivity_ocp, json_file=ocp_json_file, build=False, generate=False)
    else:
        sensitivity_solver = AcadosOcpSolver(sensitivity_ocp, json_file=ocp_json_file)

    M = chain_params_["n_mass"] - 2
    xEndRef = np.zeros((3, 1))
    xEndRef[0] = chain_params_["L"] * (M + 1) * 6
    pos0_x = np.linspace(chain_params_["xPosFirstMass"][0], xEndRef[0], M + 2)
    x0 = np.zeros((ocp.dims.nx, 1))
    x0[: 3 * (M + 1) : 3] = pos0_x[1:].reshape((M + 1, 1))

    np_test = 100

    # p_label = "L_2_0"
    # p_label = "D_2_0"
    p_label = f"C_{M}_0"

    p_idx = find_idx_for_labels(define_param_struct_symSX(chain_params_["n_mass"], disturbance=True).cat, p_label)[0]

    p_var = np.linspace(0.5 * parameter_values.cat[p_idx], 1.5 * parameter_values.cat[p_idx], np_test).flatten()

    delta_p = p_var[1] - p_var[0]

    sens_u = []
    u_opt = []

    timings_solve_ocp_solver = np.zeros((np_test))
    timings_lin_and_factorize = np.zeros((np_test))
    timings_lin_params = np.zeros((np_test))
    timings_solve_params = np.zeros((np_test))

    for i in range(np_test):
        parameter_values.cat[p_idx] = p_var[i]

        p_val = parameter_values.cat.full().flatten()
        for stage in range(ocp.dims.N + 1):
            ocp_solver.set(stage, "p", p_val)
            sensitivity_solver.set(stage, "p", p_val)

        u_opt.append(ocp_solver.solve_for_x0(x0))
        print(f"ocp_solver status {ocp_solver.status}")

        timings_solve_ocp_solver[i] = ocp_solver.get_stats("time_tot")

        ocp_solver.store_iterate(filename="iterate.json", overwrite=True, verbose=False)
        sensitivity_solver.load_iterate(filename="iterate.json", verbose=False)
        sensitivity_solver.solve_for_x0(x0, fail_on_nonzero_status=False, print_stats_on_failure=False)

        timings_lin_and_factorize[i] = sensitivity_solver.get_stats("time_tot")
        print(f"sensitivity_solver status {sensitivity_solver.status}")

        # Calculate the policy gradient
        _, sens_u_ = sensitivity_solver.eval_solution_sensitivity(0, "params_global")
        timings_lin_params[i] = sensitivity_solver.get_stats("time_solution_sens_lin")
        timings_solve_params[i] = sensitivity_solver.get_stats("time_solution_sens_solve")

        sens_u.append(sens_u_[:, p_idx])

    timing_results = {
        "NLP solve": timings_solve_ocp_solver,
        "prepare \& factorize exact Hessian QP": timings_lin_and_factorize,
        "eval rhs": timings_lin_params,
        "solve": timings_solve_params,
    }

    u_opt = np.vstack(u_opt)
    sens_u = np.vstack(sens_u)

    # Compare to numerical gradients
    sens_u_fd = np.gradient(u_opt, p_var, axis=0)
    u_opt_reconstructed_fd = np.cumsum(sens_u_fd, axis=0) * delta_p + u_opt[0, :]
    u_opt_reconstructed_fd += u_opt[0, :] - u_opt_reconstructed_fd[0, :]

    u_opt_reconstructed_acados = np.cumsum(sens_u, axis=0) * delta_p + u_opt[0, :]
    u_opt_reconstructed_acados += u_opt[0, :] - u_opt_reconstructed_acados[0, :]

    plt.figure(figsize=(7, 10))
    for col in range(3):
        plt.subplot(4, 1, col + 1)
        plt.plot(p_var, u_opt[:, col], label=f"$u^*_{col}$")
        plt.plot(p_var, u_opt_reconstructed_fd[:, col], label=f"$u^*_{col}$, reconstructed with fd gradients", linestyle="--")
        plt.plot(
            p_var, u_opt_reconstructed_acados[:, col], label=f"$u^*_{col}$, reconstructed with acados gradients", linestyle=":"
        )
        plt.ylabel(f"$u^*_{col}$")
        plt.grid(True)
        plt.legend()
        plt.xlim(p_var[0], p_var[-1])

    for col in range(3):
        plt.subplot(4, 1, 4)
        plt.plot(p_var, np.abs(sens_u[:, col] - sens_u_fd[:, col]), label=f"$u^*_{col}$", linestyle="--")

    plt.ylabel("abs difference")
    plt.grid(True)
    plt.legend()
    plt.yscale("log")
    plt.xlabel(p_label)
    plt.xlim(p_var[0], p_var[-1])

    plot_timings([timing_results], timing_results.keys(), ["acados"], figure_filename=None)

    plt.show()


if __name__ == "__main__":
    chain_params = get_chain_params()
    chain_params["n_mass"] = 3
    main_parametric(qp_solver_ric_alg=0, chain_params_=chain_params, generate_code=True)

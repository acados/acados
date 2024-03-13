"""
Test for solution sensitivities with many parameters.
"""

import os
import scipy
import numpy as np
from casadi import SX, norm_2, vertcat
import casadi as ca
from casadi.tools import struct_symSX, struct_SX, entry
from casadi.tools.structure3 import DMStruct
import matplotlib.pyplot as plt
from acados_template import AcadosModel, AcadosSim, AcadosSimSolver, AcadosOcp, AcadosOcpSolver

from typing import Tuple
from plot_utils import (
    plot_chain_control_traj,
    plot_chain_position_traj,
    plot_chain_velocity_traj,
    animate_chain_position,
)


def get_chain_params():
    """Get chain parameters."""
    params = {}

    params["n_mass"] = 5
    params["Ts"] = 0.2
    params["Tsim"] = 5
    params["N"] = 40
    params["u_init"] = np.array([-1, 1, 1])
    params["with_wall"] = True
    params["yPosWall"] = -0.05  # Dimitris: - 0.1;
    params["m"] = 0.033  # mass of the balls
    params["D"] = 1.0  # spring constant
    params["L"] = 0.033  # rest length of spring
    params["C"] = 0.1  # damping constant
    params["perturb_scale"] = 1e-2

    # params["m"] = np.random.normal(params["m_nom"], params["perturb_scale"] * params["m_nom"], params["n_mass"])
    # params["D"] = np.random.normal(params["D_nom"], params["perturb_scale"] * params["m_nom"], (params["n_mass"], 3))
    # params["L"] = np.random.normal(params["L_nom"], params["perturb_scale"] * params["m_nom"], (params["n_mass"], 3))
    # params["C"] = np.random.normal(params["C_nom"], params["perturb_scale"] * params["m_nom"], (params["n_mass"], 3))

    params["save_results"] = True
    params["show_plots"] = True
    params["nlp_iter"] = 50
    params["seed"] = 50
    params["nlp_tol"] = 1e-5

    return params


def export_chain_mass_model(n_mass: int, disturbance=False) -> Tuple[AcadosModel, DMStruct]:
    """Export chain mass model for acados."""
    x0 = np.array([0, 0, 0])  # fix mass (at wall)

    M = n_mass - 2  # number of intermediate massesu

    nx = (2 * M + 1) * 3  # differential states
    nu = 3  # control inputs

    xpos = SX.sym("xpos", (M + 1) * 3, 1)  # position of fix mass eliminated
    xvel = SX.sym("xvel", M * 3, 1)
    u = SX.sym("u", nu, 1)
    xdot = SX.sym("xdot", nx, 1)

    f = SX.zeros(3 * M, 1)  # force on intermediate masses

    n_link = n_mass - 1

    param_enries = [
        entry("m", shape=1, repeat=n_link),
        entry("D", shape=3, repeat=n_link),
        entry("L", shape=3, repeat=n_link),
        entry("C", shape=3, repeat=n_link),
    ]

    if disturbance:
        param_enries.append(entry("w", shape=3, repeat=n_mass - 2))

    p = struct_symSX(param_enries)

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
        # F = p["D", i] / p["m", i] * (1 - p["L", i] / norm_2(dist)) * dist

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
        # F = p["C", i] * vel

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

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p.cat
    model.name = model_name

    p_map = p(0)

    return model, p_map


def compute_steady_state(model, p, xPosFirstMass, xEndRef):
    """Compute steady state for chain mass model."""

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

    wrest = sol["x"].full()
    xrest = wrest[:nx]

    return xrest


def sampleFromEllipsoid(w, Z):
    """
    draws uniform sample from ellipsoid with center w and variability matrix Z
    """

    n = w.shape[0]  # dimension
    lam, v = np.linalg.eig(Z)

    # sample in hypersphere
    r = np.random.rand() ** (1 / n)  # radial position of sample
    x = np.random.randn(n)
    x = x / np.linalg.norm(x)
    x *= r
    # project to ellipsoid
    y = v @ (np.sqrt(lam) * x) + w

    return y


def export_chain_mass_integrator(n_mass, p_, disturbance=True) -> AcadosSimSolver:
    """Export chain mass integrator for acados."""
    sim = AcadosSim()
    # simulation options
    Ts = 0.2

    # export model
    model, _ = export_chain_mass_model(n_mass, disturbance=disturbance)

    # set model
    sim.model = model

    # disturbances
    # nparam = 3 * M
    sim.parameter_values = p_.cat.full().flatten()

    # solver options
    sim.solver_options.integrator_type = "IRK"

    sim.solver_options.num_stages = 2
    sim.solver_options.num_steps = 2
    # sim.solver_options.nlp_solver_tol_eq = 1e-9

    # set prediction horizon
    sim.solver_options.T = Ts

    # acados_ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp_' + model.name + '.json')
    acados_integrator = AcadosSimSolver(sim, json_file="acados_ocp_" + model.name + ".json")

    return acados_integrator


def main(_chain_params: dict):
    """Main function."""
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # chain parameters
    n_mass = _chain_params["n_mass"]
    M = _chain_params["n_mass"] - 2  # number of intermediate masses
    Ts = _chain_params["Ts"]
    Tsim = _chain_params["Tsim"]
    N = _chain_params["N"]
    u_init = _chain_params["u_init"]
    with_wall = _chain_params["with_wall"]
    yPosWall = _chain_params["yPosWall"]

    perturb_scale = _chain_params["perturb_scale"]
    W = perturb_scale * np.eye(M)

    m = [np.random.normal(_chain_params["m"], 0.1 * _chain_params["m"]) for _ in range(n_mass - 1)]
    D = [np.random.normal(_chain_params["D"], 0.1 * _chain_params["D"], 3) for _ in range(n_mass - 1)]
    L = [np.random.normal(_chain_params["L"], 0.1 * _chain_params["L"], 3) for _ in range(n_mass - 1)]
    C = [np.random.normal(_chain_params["C"], 0.1 * _chain_params["C"], 3) for _ in range(n_mass - 1)]

    nlp_iter = _chain_params["nlp_iter"]
    nlp_tol = _chain_params["nlp_tol"]
    save_results = _chain_params["save_results"]
    show_plots = _chain_params["show_plots"]
    seed = _chain_params["seed"]

    np.random.seed(seed)

    # nparam = 3 * M
    # w = perturb_scale * np.eye(nparam)

    # export model
    model, p = export_chain_mass_model(n_mass, disturbance=True)

    # Inermediate masses
    for i_mass in range(n_mass - 1):
        p["m", i_mass] = m[i_mass]
    for i_mass in range(n_mass - 2):
        p["w", i_mass] = np.zeros(3)

    # Links
    for i_link in range(n_mass - 1):
        p["D", i_link] = D[i_link]
        p["L", i_link] = L[i_link]
        p["C", i_link] = C[i_link]

    # set model
    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx
    Tf = N * Ts

    # initial state
    xPosFirstMass = np.zeros((3, 1))

    xEndRef = np.zeros((3, 1))
    xEndRef[0] = _chain_params["L"] * (M + 1) * 6
    xEndRef[2] = 1.0
    # xEndRef[0] = M * chain_params["L"]

    # xrest = compute_steady_state(n_mass, p, xPosFirstMass, xEndRef, disturbance=True)
    xrest = compute_steady_state(model, p, xPosFirstMass, xEndRef)

    x0 = xrest

    if True:
        pos = np.zeros((n_mass, 3))
        pos[0, :] = xPosFirstMass.T
        for i in range(1, n_mass):
            j = i - 1
            pos[i, :] = x0[j * 3 : (j + 1) * 3].T

        plt.figure()
        plt.plot(pos[:, 0], pos[:, 2], "o-")
        plt.xlabel("x")
        plt.ylabel("z")
        plt.title("Initial state")
        plt.grid(True)
        plt.show()

        print(f"x0 = {x0}")

    # set dimensions
    ocp.dims.N = N

    # set cost module
    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"

    Q = 2 * np.diagflat(np.ones((nx, 1)))
    q_diag = np.ones((nx, 1))
    strong_penalty = M + 1
    q_diag[3 * M] = strong_penalty
    q_diag[3 * M + 1] = strong_penalty
    q_diag[3 * M + 2] = strong_penalty
    Q = 2 * np.diagflat(q_diag)

    R = 2 * np.diagflat(1e-2 * np.ones((nu, 1)))

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Q

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx, :nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[nx : nx + nu, :] = np.eye(nu)
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    yref = np.vstack((xrest, np.zeros((nu, 1)))).flatten()
    ocp.cost.yref = yref
    ocp.cost.yref_e = xrest.flatten()

    # set constraints
    umax = 1 * np.ones((nu,))

    ocp.constraints.lbu = -umax
    ocp.constraints.ubu = umax
    ocp.constraints.x0 = x0.reshape((nx,))
    ocp.constraints.idxbu = np.array(range(nu))

    # disturbances
    ocp.parameter_values = p.cat.full().flatten()

    # wall constraint
    if with_wall:
        nbx = M + 1
        Jbx = np.zeros((nbx, nx))
        for i in range(nbx):
            Jbx[i, 3 * i + 1] = 1.0

        ocp.constraints.Jbx = Jbx
        ocp.constraints.lbx = yPosWall * np.ones((nbx,))
        ocp.constraints.ubx = 1e9 * np.ones((nbx,))

        # slacks
        ocp.constraints.Jsbx = np.eye(nbx)
        L2_pen = 1e3
        L1_pen = 1
        ocp.cost.Zl = L2_pen * np.ones((nbx,))
        ocp.cost.Zu = L2_pen * np.ones((nbx,))
        ocp.cost.zl = L1_pen * np.ones((nbx,))
        ocp.cost.zu = L1_pen * np.ones((nbx,))

    # ERK4 integrator

    # solver options
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"  # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "IRK"
    ocp.solver_options.nlp_solver_type = "SQP"  # SQP_RTI
    ocp.solver_options.nlp_solver_max_iter = nlp_iter

    ocp.solver_options.sim_method_num_stages = 2
    ocp.solver_options.sim_method_num_steps = 2
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.qp_tol = nlp_tol
    ocp.solver_options.tol = nlp_tol
    # ocp.solver_options.nlp_solver_tol_eq = 1e-9

    # set prediction horizon
    ocp.solver_options.tf = Tf

    acados_ocp_solver = AcadosOcpSolver(ocp, json_file="acados_ocp_" + model.name + ".json")

    # acados_integrator = AcadosSimSolver(ocp, json_file = 'acados_ocp_' + model.name + '.json')
    acados_integrator = export_chain_mass_integrator(n_mass, p)

    # %% get initial state from xrest
    xcurrent = x0.reshape((nx,))
    for i in range(5):
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", u_init)

        status = acados_integrator.solve()
        if status != 0:
            raise Exception("acados integrator returned status {}. Exiting.".format(status))

        # update state
        xcurrent = acados_integrator.get("x")

    print(f"Initial state: {xcurrent}")

    # %% actual simulation
    N_sim = int(np.floor(Tsim / Ts))
    simX = np.ndarray((N_sim + 1, nx))
    simU = np.ndarray((N_sim, nu))
    wall_dist = np.zeros((N_sim,))

    timings = np.zeros((N_sim,))

    simX[0, :] = xcurrent

    # closed loop
    for i in range(N_sim):
        # solve ocp
        acados_ocp_solver.set(0, "lbx", xcurrent)
        acados_ocp_solver.set(0, "ubx", xcurrent)

        status = acados_ocp_solver.solve()
        timings[i] = acados_ocp_solver.get_stats("time_tot")

        if status != 0:
            raise Exception("acados acados_ocp_solver returned status {} in time step {}. Exiting.".format(status, i))

        simU[i, :] = acados_ocp_solver.get(0, "u")
        print("control at time", i, ":", simU[i, :])

        # simulate system
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", simU[i, :])

        p["w"] = sampleFromEllipsoid(np.zeros((W.shape[0],)), W)

        acados_integrator.set("p", p.cat.full().flatten())

        status = acados_integrator.solve()
        if status != 0:
            raise Exception("acados integrator returned status {}. Exiting.".format(status))

        # update state
        xcurrent = acados_integrator.get("x")
        simX[i + 1, :] = xcurrent

        # xOcpPredict = acados_ocp_solver.get(1, "x")
        # print("model mismatch = ", str(np.max(xOcpPredict - xcurrent)))
        yPos = xcurrent[range(1, 3 * M + 1, 3)]
        wall_dist[i] = np.min(yPos - yPosWall)
        print("time i = ", str(i), " dist2wall ", str(wall_dist[i]))

    print("dist2wall (minimum over simulation) ", str(np.min(wall_dist)))

    if os.environ.get("ACADOS_ON_CI") is None and show_plots:
        plot_chain_control_traj(simU)
        plot_chain_position_traj(simX, yPosWall=yPosWall)
        plot_chain_velocity_traj(simX)

        animate_chain_position(simX, xPosFirstMass, yPosWall=yPosWall)
        # animate_chain_position_3D(simX, xPosFirstMass)

        plt.show()


if __name__ == "__main__":
    chain_params = get_chain_params()

    main(chain_params)

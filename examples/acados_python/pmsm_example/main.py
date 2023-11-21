from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver
from casadi import SX, vertcat, sqrt, cos, sin
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.linalg


FORMULATION = 1
# Tracking MPC
# ===================
# 0 for hexagon (no terminalset)

# Economic MPC
# ===================
# 1 PDC Voltage + PDC Terminal Set

i_d_ref = -125  # Setpoints only valid for Formulation 0
i_q_ref = 10
w_val = 2000.0  # do not change in this script, some of
# the values below are calculate with a fix w_val = 2000 1/S
tau_wal = 10.0


# constants
L_d = 107e-6
L_q = 150e-6
R_m = 18.15e-3
K_m = 13.8e-3
N_P = 5.0
u_max = 48
i_max = 155.0  # only for simulation

phi = 0.0

x0Start = np.array([0, 0])

N = 2
Ts_sim = 125e-6

WARMSTART_ITERS = 200

INPUT_REG = 1e-2

# setpoint MPC with hexagon (current slack)
if FORMULATION == 0:  # (works)
    Weight_TUNING = 1
    Weight_E_TUNING = 1
    SLACK_TUNING = 0
    SLACK_E_TUNING = 0
    wd = 1
    wq = 1
    Ts = 0.000125
    SLACK_TUNINGHessian = 0

# full EMPC
if FORMULATION == 1:  # (works)
    Weight_TUNING = 1e-1
    Weight_E_TUNING = 1e-1
    SLACK_TUNING = 1e3
    SLACK_E_TUNING = 1e3
    wd = 1
    wq = 1
    Ts = 0.000250
    SLACK_TUNINGHessian = 0

# Model of the PMSM
# ====================================================================
def export_pmsm_model():

    model_name = "pmsm"

    # set up states
    i_d = SX.sym("i_d")
    i_q = SX.sym("i_q")
    x = vertcat(i_d, i_q)

    # set up controls
    u_d = SX.sym("u_d")
    u_q = SX.sym("u_q")
    u = vertcat(u_d, u_q)

    # set up xdot
    i_d_dot = SX.sym("i_d_dot")
    i_q_dot = SX.sym("i_q_dot")
    xdot = vertcat(i_d_dot, i_q_dot)

    # set up parameters
    omega = SX.sym("omega")  # speed
    dist_d = SX.sym("dist_d")  # d disturbance
    dist_q = SX.sym("dist_q")  # q disturbance
    tau_des = SX.sym("tau_des")  # Desired tau
    p = vertcat(omega, dist_d, dist_q, tau_des)

    # dynamics
    fexpl = vertcat(
        -(R_m / L_d) * i_d + (L_q / L_d) * omega * i_q + u_d / L_d,
        -(L_d / L_q) * omega * i_d
        - (R_m / L_q) * i_q
        + u_q / L_q
        - (omega * K_m) / L_q,
    )

    fimpl = vertcat(
        -i_d_dot - (R_m / L_d) * i_d + (L_q / L_d) * omega * i_q + u_d / L_d,
        -i_q_dot
        - (L_d / L_q) * omega * i_d
        - (R_m / L_q) * i_q
        + u_q / L_q
        - (omega * K_m) / L_q,
    )

    model = AcadosModel()

    if FORMULATION == 1:
        # export_torqueline_pd():
        # torque and voltage constraints
        r = SX.sym("r", 3, 1)
        model.con_r_in_phi = r
        model.con_phi_expr = vertcat(r[0], r[1] ** 2 + r[2] ** 2)
        model.con_r_expr = vertcat(
            tau_des - 1.5 * N_P * ((L_d - L_q) * i_d * i_q + K_m * i_q), u_d, u_q
        )

        # also at initial node
        model.con_r_in_phi_0 = model.con_r_in_phi
        model.con_phi_expr_0 = model.con_phi_expr
        model.con_r_expr_0 = model.con_r_expr

        # export_torquelineEnd_pd():
        alpha = R_m**2 + omega**2 * L_d**2
        beta = R_m**2 + omega**2 * L_q**2
        gamma = 2 * R_m * omega * (L_d - L_q)
        delta = 2 * omega**2 * L_d * K_m
        epsilon = 2 * R_m * omega * K_m
        rho = omega**2 * K_m**2

        # torque and voltage constraints
        r_e = SX.sym("r_e", 3, 1)
        model.con_phi_expr_e = vertcat(
            r_e[0],
            alpha * r_e[1] ** 2
            + beta * r_e[2] ** 2
            + gamma * r_e[1] * r_e[2]
            + delta * r_e[1]
            + epsilon * r_e[2]
            + rho,
        )
        model.con_r_expr_e = vertcat(
            tau_des - 1.5 * N_P * ((L_d - L_q) * i_d * i_q + K_m * i_q), i_d, i_q
        )
        model.con_r_in_phi_e = r_e

    model.f_impl_expr = fimpl
    model.f_expl_expr = fexpl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = []
    model.p = p
    model.name = model_name

    return model


# Hexagon Voltage Constraints
# ====================================================================
def get_general_constraints_DC():

    # polytopic constraint on the input
    s3 = sqrt(3)
    cphi = cos(phi)
    sphi = sin(phi)

    h1 = s3 * cphi + sphi
    h2 = sphi
    h3 = -s3 * cphi + sphi
    h7 = -s3 * sphi + cphi
    h8 = cphi
    h9 = s3 * sphi + cphi

    # form D and C matrices
    # (acados C interface works with column major format)
    D = np.array([[h1, h7], [h2, h8], [h3, h9]])
    C = np.array([[0, 0], [0, 0], [0, 0]])

    g1 = 2.0 / s3 * u_max
    g2 = 1.0 / s3 * u_max

    lg = np.array([-g1, -g2, -g1])
    ug = np.array([g1, g2, g1])

    res = dict()
    res["D"] = D
    res["C"] = C
    res["lg"] = lg
    res["ug"] = ug

    return res


# Hexagon Terminal Constraints
# ====================================================================
def get_general_terminal_constraints_DC():

    # polytopic constraint on the input
    s3 = sqrt(3)
    cphi = cos(phi)
    sphi = sin(phi)

    h1 = s3 * cphi + sphi
    h2 = sphi
    h3 = -s3 * cphi + sphi
    h7 = -s3 * sphi + cphi
    h8 = cphi
    h9 = s3 * sphi + cphi

    # form D and C matrices
    D = np.array([[h1, h7], [h2, h8], [h3, h9]])

    A = np.array([[-R_m / L_d, w_val * L_q / L_d], [-w_val * L_d / L_q, -R_m / L_q]])
    invB = np.array([[L_d, 0], [0, L_q]])
    f = np.array([[0], [-K_m * w_val / L_q]])

    invBA = invB.dot(A)
    Ce = D.dot(invBA)

    g1 = 2.0 / s3 * u_max
    g2 = 1.0 / s3 * u_max

    invBf = invB.dot(f)
    lge = -np.array([[g1], [g2], [g1]]) - D.dot(invBf)
    uge = +np.array([[g1], [g2], [g1]]) - D.dot(invBf)

    res = dict()
    res["Ce"] = Ce
    res["lge"] = lge.flatten()
    res["uge"] = uge.flatten()

    return res


def run_simulation(qp_solver="FULL_CONDENSING_HPIPM", show_plots=False, verbose=False):
    # create render arguments
    ocp = AcadosOcp()

    # export model
    model = export_pmsm_model()
    ocp.model = model

    # model dims
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nu + nx
    ny_e = nx
    Tf = N * Ts

    # set ocp_nlp_dimensions
    ocp.dims.N = N
    # NOTE: all dimensions but N are now detected by the interface.

    # set weighting matrices
    nlp_cost = ocp.cost

    Q = np.eye(nx)
    Q[0, 0] = wd * Weight_TUNING
    Q[1, 1] = wq * Weight_TUNING

    Q_e = np.eye(nx)
    Q_e[0, 0] = wd * Weight_E_TUNING * Tf / N
    Q_e[1, 1] = wq * Weight_E_TUNING * Tf / N

    R = np.eye(nu)
    R[0, 0] = INPUT_REG
    R[1, 1] = INPUT_REG

    nlp_cost.W = scipy.linalg.block_diag(Q, R)  # weight matrix
    nlp_cost.W_e = Q_e  # weight matrix for Mayer term

    Vu = np.zeros((ny, nu))
    Vu[2, 0] = 1.0
    Vu[3, 1] = 1.0
    nlp_cost.Vu = Vu

    Vx = np.zeros((ny, nx))
    Vx[0, 0] = 1.0
    Vx[1, 1] = 1.0
    nlp_cost.Vx = Vx

    Vx_e = np.zeros((ny_e, nx))
    Vx_e[0, 0] = 1.0
    Vx_e[1, 1] = 1.0
    nlp_cost.Vx_e = Vx_e

    if FORMULATION == 0:
        nlp_cost.yref = np.zeros((ny,))
        nlp_cost.yref[0] = i_d_ref
        nlp_cost.yref[1] = i_q_ref
        nlp_cost.yref[2] = 0
        nlp_cost.yref[3] = 0
        nlp_cost.yref_e = np.zeros((ny_e,))
        nlp_cost.yref_e[0] = i_d_ref
        nlp_cost.yref_e[1] = i_q_ref
    else:
        nlp_cost.yref = np.zeros((ny,))
        nlp_cost.yref_e = np.zeros((ny_e,))

    # get D and C
    res = get_general_constraints_DC()
    D = res["D"]
    C = res["C"]
    lg = res["lg"]
    ug = res["ug"]

    res = get_general_terminal_constraints_DC()
    Ce = res["Ce"]
    lge = res["lge"]
    uge = res["uge"]

    # setting bounds
    nlp_con = ocp.constraints

    if FORMULATION == 0:
        # setting general constraints --> lg <= D*u + C*u <= ug
        nlp_con.D = D
        nlp_con.C = C
        nlp_con.lg = lg
        nlp_con.ug = ug

        nlp_con.x0 = x0Start

    if FORMULATION == 1:
        # setting general constraints --> lg <= D*u + C*u <= ug
        nlp_con.D = D
        nlp_con.C = C
        nlp_con.lg = lg
        nlp_con.ug = ug
        nlp_con.C_e = Ce
        nlp_con.lg_e = lge
        nlp_con.ug_e = uge

        # lower gradient/hessian wrt lower slack
        nlp_cost.zl = SLACK_TUNING * np.array([1])
        nlp_cost.Zl = SLACK_TUNINGHessian * np.ones((1,))  # hessian
        nlp_cost.zu = SLACK_TUNING * np.array([1])
        nlp_cost.Zu = SLACK_TUNINGHessian * np.ones((1,))  # hessian
        # _e
        nlp_cost.zl_e = SLACK_E_TUNING * np.array([1])
        nlp_cost.Zl_e = SLACK_TUNINGHessian * np.ones((1,))  # hessian
        nlp_cost.zu_e = SLACK_E_TUNING * np.array([1])
        nlp_cost.Zu_e = SLACK_TUNINGHessian * np.ones((1,))  # hessian
        # _0
        nlp_cost.zl_0 = SLACK_E_TUNING * np.array([1])
        nlp_cost.Zl_0 = SLACK_TUNINGHessian * np.ones((1,))  # hessian
        nlp_cost.zu_0 = SLACK_E_TUNING * np.array([1])
        nlp_cost.Zu_0 = SLACK_TUNINGHessian * np.ones((1,))  # hessian


        lphi = np.array([0, -1e9]
        )  # 1st torque constraint | 2nd voltage constraint
        uphi = np.array([0, u_max**2 / 3])

        nlp_con.lphi = lphi
        nlp_con.lphi_0 = lphi
        nlp_con.lphi_e = lphi

        nlp_con.uphi = uphi
        nlp_con.uphi_0 = uphi
        nlp_con.uphi_e = uphi

        nlp_con.idxsphi = np.array([0])
        nlp_con.idxsphi_e = np.array([0])
        nlp_con.idxsphi_0 = np.array([0])

        nlp_con.x0 = x0Start

    # setting parameters
    ocp.parameter_values = np.array([w_val, 0.0, 0.0, tau_wal])

    # set QP solver
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "IRK"
    # ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.sim_method_num_stages = 1

    # ocp.solver_options.qp_solver_tol_stat = 1e-4
    ocp.solver_options.qp_solver_tol_eq = 1e-4
    ocp.solver_options.qp_solver_tol_ineq = 1e-4
    ocp.solver_options.qp_solver_tol_comp = 1e-4

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = "SQP_RTI"
    # ocp.solver_options.nlp_solver_type = 'SQP'

    file_name = "acados_ocp.json"

    acados_solver = AcadosOcpSolver(ocp, json_file=file_name)

    # test setter:
    for i in range(1,N):
        acados_solver.cost_set(i, 'zl', ocp.cost.zl)
        acados_solver.cost_set(i, 'zu', ocp.cost.zu)
        acados_solver.cost_set(i, 'Zl', ocp.cost.Zl)
        acados_solver.cost_set(i, 'Zu', ocp.cost.Zu)
        acados_solver.cost_set(i, 'W', ocp.cost.W)
        acados_solver.cost_set(i, 'yref', ocp.cost.yref)

    # test constraints set with polytopic constraints
    acados_solver.constraints_set(1, "lg", lg)
    acados_solver.constraints_set(1, "ug", ug)
    acados_solver.constraints_set(1, "D", D, api="new")
    acados_solver.constraints_set(1, "C", C, api="new")

    # closed loop simulation
    Nsim = 20

    simX = np.ndarray((Nsim, nx))
    simU = np.ndarray((Nsim, nu))
    simXR = np.ndarray((Nsim + 1, nx))
    simXRN = np.ndarray((Nsim, nx))

    print("=========================================")
    print("Mode = ", FORMULATION)
    print("speed_el = ", w_val)
    print("Sample Time = ", Ts)
    print("=========================================")
    print("Initial Condition")
    print("id0 = ", x0Start[0])
    print("iq0 = ", x0Start[1])

    # get initial condition for real simulation
    simXR[0, 0] = x0Start[0]
    simXR[0, 1] = x0Start[1]

    xvec = np.array([[x0Start[0]], [x0Start[1]]])

    # compute warm-start
    for i in range(WARMSTART_ITERS):
        status = acados_solver.solve()

    for i in range(Nsim):
        if verbose:
            print("\n")
            print("SimulationStep = ", i)
            print("=================")

        # set options
        acados_solver.options_set("print_level", 0)
        status = acados_solver.solve()

        if status != 0:
            acados_solver.print_statistics()
            raise Exception(f"acados returned status {status} in simulation step {i}.")

        # get solution
        x0 = acados_solver.get(0, "x")
        xN = acados_solver.get(N, "x")
        u0 = acados_solver.get(0, "u")

        CPU_time = acados_solver.get_stats("time_tot")

        for j in range(nx):
            simX[i, j] = x0[j]
            simXRN[i, j] = xN[j]

        for j in range(nu):
            simU[i, j] = u0[j]

        uvec = np.array([[u0[0]], [u0[1]]])

        # real Simulation
        A = np.array(
            [[-R_m / L_d, w_val * L_q / L_d], [-w_val * L_d / L_q, -R_m / L_q]],
            dtype=float,
        )
        B = np.array([[1.0 / L_d, 0.0], [0.0, 1.0 / L_q]], dtype=float)
        f = np.array([[0.0], [-K_m * w_val / L_q]], dtype=float)

        # Euler
        # ========================
        # xvec = xvec + Ts*A*xvec + Ts*B*uvec + Ts*f

        # Heun
        # ========================
        # xvec = xvec + Ts*A*xvec + Ts*B*uvec + Ts*f + 0.5*Ts*Ts*A*A*xvec + 0.5*Ts*Ts*A*B*uvec + 0.5*Ts*Ts*A*f

        # Z-Transformation
        # ========================
        Ad = scipy.linalg.expm(A * Ts_sim)
        invA = np.linalg.inv(A)
        Bd = np.dot(invA, np.dot((Ad - np.eye(2)), B))
        fd = np.dot(invA, np.dot((Ad - np.eye(2)), f))
        xvec = np.dot(Ad, xvec) + np.dot(Bd, uvec) + fd
        # xvec_arg = np.zeros((2,0))
        xvec_arg = np.zeros((2,))
        xvec_arg[0] = xvec[0, 0]
        xvec_arg[1] = xvec[1, 0]

        if verbose:
            print("States= ", xvec)
            print("Controls= ", uvec)
            print("\n")

        # update initial condition xk+1
        acados_solver.constraints_set(0, "lbx", xvec_arg)
        acados_solver.constraints_set(0, "ubx", xvec_arg)

        for j in range(N):
            acados_solver.cost_set(j, "W", nlp_cost.W)
        acados_solver.cost_set(N, "W", nlp_cost.W_e)

        simXR[i + 1, 0] = xvec[0]
        simXR[i + 1, 1] = xvec[1]

    del acados_solver

    # tests
    xref_sol = np.array([[-77.39723565], [77.84933005]])
    err_x_final = np.max(np.abs(xvec - xref_sol))
    print(f"error in x wrt reference solution: {err_x_final:.2e}")

    uref_sol = np.array([[-24.75529731], [ 12.44582214]])
    err_u_final = np.max(np.abs(uvec - uref_sol))
    print(f"error in u wrt reference solution: {err_u_final:.2e}")

    if err_x_final > 1e-3:
        raise Exception("error wrt reference solution should be < 1e-3.")
    if err_u_final > 1e-3:
        raise Exception("error wrt reference solution should be < 1e-3.")

    # avoid plotting when running on Travis
    if os.environ.get("ACADOS_ON_CI") is None and show_plots:
        # plot results
        t = np.linspace(0.0, Ts * Nsim, Nsim)
        plt.subplot(4, 1, 1)
        plt.step(t, simU[:, 0], "r")
        plt.step(t, simU[:, 0], "ro")
        plt.plot([0, Ts * Nsim], [nlp_cost.yref[2], nlp_cost.yref[2]], "--")
        plt.title("closed-loop simulation")
        plt.ylabel("u_d")
        plt.xlabel("t")
        plt.grid(True)
        plt.subplot(4, 1, 2)
        plt.step(t, simU[:, 1], "r")
        plt.step(t, simU[:, 1], "ro")
        plt.plot([0, Ts * Nsim], [nlp_cost.yref[3], nlp_cost.yref[3]], "--")
        plt.ylabel("u_q")
        plt.xlabel("t")
        plt.grid(True)
        plt.subplot(4, 1, 3)
        plt.plot(t, simX[:, 0])
        plt.ylabel("x_d")
        plt.xlabel("t")
        plt.grid(True)
        plt.subplot(4, 1, 4)
        plt.plot(t, simX[:, 1])
        plt.ylabel("x_q")
        plt.xlabel("t")
        plt.grid(True)

        # plot hexagon
        r = 2 / 3 * u_max
        x1 = r
        y1 = 0
        x2 = r * cos(pi / 3)
        y2 = r * sin(pi / 3)
        q1 = -(y2 - y1 / x1 * x2) / (1 - x2 / x1)
        m1 = -(y1 + q1) / x1

        # box constraints
        q2 = r * sin(pi / 3)
        # -q2 <= uq  <= q2

        plt.figure()
        plt.plot(simU[:, 0], simU[:, 1], "o")
        plt.xlabel("ud")
        plt.ylabel("uq")
        ud = np.linspace(-1.5 * u_max, 1.5 * u_max, 100)
        plt.plot(ud, -m1 * ud - q1)
        plt.plot(ud, -m1 * ud + q1)
        plt.plot(ud, +m1 * ud - q1)
        plt.plot(ud, +m1 * ud + q1)
        plt.plot(ud, -q2 * np.ones((100, 1)))
        plt.plot(ud, q2 * np.ones((100, 1)))
        plt.grid(True)
        ax = plt.gca()
        ax.set_xlim([-u_max, u_max])
        ax.set_ylim([-u_max, u_max])
        circle = plt.Circle((0, 0), u_max / np.sqrt(3), color="red", fill=False)
        ax.add_artist(circle)

        delta = 100
        x = np.linspace(-i_max, i_max / 3, delta)
        y = np.linspace(-i_max, i_max, delta)
        XV, YV = np.meshgrid(x, y)

        alpha = R_m**2 + w_val**2 * L_d**2
        beta = R_m**2 + w_val**2 * L_q**2
        gamma = 2 * R_m * w_val * (L_d - L_q)
        delta = 2 * w_val**2 * L_d * K_m
        epsilon = 2 * R_m * w_val * K_m
        rho = w_val**2 * K_m**2 - (u_max) ** 2 / 3

        FeaSetV = (
            alpha * XV**2
            + beta * YV**2
            + gamma * XV * YV
            + delta * XV
            + epsilon * YV
            + rho
        )
        TauV = 1.5 * N_P * ((L_d - L_q) * XV * YV + K_m * YV)

        # trajectory in the dq-plane
        plt.figure()
        plt.plot(simXR[:, 0], simXR[:, 1], "bo")
        plt.plot(simXR[:, 0], simXR[:, 1], "b")
        plt.plot(simX[:, 0], simX[:, 1], "r")
        plt.plot(simX[:, 0], simX[:, 1], "r*")
        plt.contour(XV, YV, FeaSetV, [0], colors="r")
        plt.contour(XV, YV, TauV, [tau_wal], colors="g")
        plt.xlabel("x_d")
        plt.ylabel("x_q")
        plt.grid(True)

        S1 = -0.24543672 * XV + 0.50146524 * YV
        S2 = -0.214 * XV - 0.01815 * YV
        S3 = -0.18256328 * XV - 0.53776524 * YV

        plt.figure()
        plt.plot(simXRN[:, 0], simXRN[:, 1], "bo")
        plt.plot(simXRN[:, 0], simXRN[:, 1], "b")
        plt.contour(XV, YV, FeaSetV, [0], colors="r")
        plt.contour(XV, YV, S1, [-27.82562584, 83.02562584], colors="k")
        plt.contour(XV, YV, S2, [-0.11281292, 55.31281292], colors="k")
        plt.contour(XV, YV, S3, [-27.82562584, 83.02562584], colors="k")
        plt.contour(XV, YV, TauV, [tau_wal], colors="g")
        plt.xlabel("x_d")
        plt.ylabel("x_q")
        plt.grid(True)

        plt.show()


def main():
    for qp_solver in ["FULL_CONDENSING_DAQP", "FULL_CONDENSING_HPIPM"]: # , "FULL_CONDENSING_QPOASES"
        run_simulation(qp_solver=qp_solver)


if __name__ == "__main__":
    main()

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver
from robot_model import export_robot_model
import numpy as np
import scipy.linalg
from utils import plot_robot


X0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])  # Intitalize the states [x,y,v,th,th_d]
N_horizon = 50  # Define the number of discretization steps
T_horizon = 2.0  # Define the prediction horizon
F_max = 10  # Define the max force allowed
F_min = -10
T_max = 1.0
T_max = -1.0


def create_ocp_solver_description() -> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    model = export_robot_model()
    ocp.model = model
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu

    # set dimensions
    ocp.dims.N = N_horizon

    # set cost
    Q_mat = 2 * np.diag([1e3, 1e3, 0.0, 0.0, 0.0])  # [x,y,x_d,y_d,th,th_d]
    R_mat = 2 * 5 * np.diag([1e-1, 1e-2])

    ocp.cost.cost_type = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"

    ny = nx + nu
    ny_e = nx

    ocp.cost.W_e = Q_mat
    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx, :nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[nx : nx + nu, 0:nu] = np.eye(nu)
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    # set constraints
    ocp.constraints.lbu = np.array([-F_max])
    ocp.constraints.ubu = np.array([+F_max])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = X0

    # set options
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"  # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"  # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = "IRK"
    ocp.solver_options.nlp_solver_type = "SQP"  # SQP_RTI, SQP
    ocp.solver_options.nlp_solver_max_iter = 400
    # ocp.solver_options.levenberg_marquardt = 1e-2

    # set prediction horizon
    ocp.solver_options.tf = T_horizon

    return ocp


def closed_loop_simulation():

    # create solvers
    ocp = create_ocp_solver_description()
    acados_ocp_solver = AcadosOcpSolver(
        ocp, json_file="acados_ocp_" + ocp.model.name + ".json"
    )
    acados_integrator = AcadosSimSolver(
        ocp, json_file="acados_ocp_" + ocp.model.name + ".json"
    )

    # prepare simulation
    Nsim = 100
    nx = ocp.model.x.size()[0]
    nu = ocp.model.u.size()[0]

    simX = np.ndarray((Nsim + 1, nx))
    simU = np.ndarray((Nsim, nu))

    xcurrent = X0
    simX[0, :] = xcurrent

    # initialize solver
    for stage in range(N_horizon + 1):
        acados_ocp_solver.set(stage, "x", 0.0 * np.ones(xcurrent.shape))
    for stage in range(N_horizon):
        acados_ocp_solver.set(stage, "u", np.zeros((nu,)))

    # closed loop
    for i in range(Nsim):

        # set initial state constraint
        acados_ocp_solver.set(0, "lbx", xcurrent)
        acados_ocp_solver.set(0, "ubx", xcurrent)

        # update yref
        for j in range(N_horizon):
            yref = np.array([1, 1, 0, 0, 0, 0, 0])
            acados_ocp_solver.set(j, "yref", yref)
        yref_N = np.array([1, 1, 0, 0, 0])
        acados_ocp_solver.set(N_horizon, "yref", yref_N)

        # solve ocp
        status = acados_ocp_solver.solve()

        if status not in [0, 2]:
            acados_ocp_solver.print_statistics()
            plot_robot(
                np.linspace(0, T_horizon / N_horizon * i, i + 1),
                F_max,
                simU[:i, :],
                simX[: i + 1, :],
            )
            raise Exception(
                f"acados acados_ocp_solver returned status {status} in closed loop instance {i} with {xcurrent}"
            )

        if status == 2:
            print(
                f"acados acados_ocp_solver returned status {status} in closed loop instance {i} with {xcurrent}"
            )
        simU[i, :] = acados_ocp_solver.get(0, "u")

        # simulate system
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", simU[i, :])

        status = acados_integrator.solve()
        if status != 0:
            raise Exception(
                f"acados integrator returned status {status} in closed loop instance {i}"
            )

        # update state
        xcurrent = acados_integrator.get("x")
        simX[i + 1, :] = xcurrent

    # plot results
    plot_robot(
        np.linspace(0, T_horizon / N_horizon * Nsim, Nsim + 1), [F_max, None], simU, simX
    )


if __name__ == "__main__":
    closed_loop_simulation()

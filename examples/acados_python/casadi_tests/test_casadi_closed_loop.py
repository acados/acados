import sys
sys.path.insert(0, '../getting_started')

import numpy as np
import casadi as ca

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver, AcadosSim, AcadosSimSolver
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum

PLOT = False

def formulate_ocp_and_sim(Tf: float = 1.0, N: int = 20, dt: float = 0.05)-> tuple[AcadosOcp, AcadosSimSolver]:
    ### create ocp
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    
    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0, np.pi, 0, 0])  # initial state
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization

    ### set integrator
    sim = AcadosSim()
    sim.model = model
    sim.solver_options.T = dt
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 2
    sim.solver_options.integrator_type = 'ERK'

    integrator = AcadosSimSolver(sim, verbose=False)

    return ocp, integrator

def main():
    Tsim = 2.0
    Tf = 0.5
    dt = 0.05
    Nsim = int(Tsim/dt)
    N_horizon = int(Tf/dt)
    ocp, integrator = formulate_ocp_and_sim(Tf, N_horizon, dt)
    ocp.make_consistent()

    simX = np.zeros((Nsim+1, ocp.dims.nx))
    simU = np.zeros((Nsim, ocp.dims.nu))
    simX_casadi = np.zeros((Nsim+1, ocp.dims.nx))
    simU_casadi = np.zeros((Nsim, ocp.dims.nu))
    simX[0,:] = ocp.constraints.x0
    simX_casadi[0,:] = ocp.constraints.x0

    # steady-state
    xs = np.array([[0.0, 0.0, 0.0, 0.0]]).T
    us = np.array([[0.0]]).T

    # constant ref
    X_ref = np.tile(xs, Nsim + 1).T
    U_ref = np.tile(us, Nsim).T

    # reference jump
    xs2 = np.array([0.5, 0, 0, 0])
    Njump = int(Nsim / 4)
    X_ref[Njump : 3*Njump, :] = xs2

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp,verbose=False)
    # do some initial iterations to start with a good initial guess
    num_iter_initial = 10
    for _ in range(num_iter_initial):
        ocp_solver.solve_for_x0(x0_bar = simX[0,:], fail_on_nonzero_status=False)
    # closed loop with acados
    for i in range(Nsim):
        yref = np.concatenate((X_ref[i, :], U_ref[i, :]))
        for stage in range(ocp.solver_options.N_horizon):
            ocp_solver.set(stage, "yref", yref)
        ocp_solver.set(ocp.solver_options.N_horizon, "yref", X_ref[i, :])
        simU[i, :] = ocp_solver.solve_for_x0(simX[i, :])
        simX[i+1, :] = integrator.simulate(simX[i, :], simU[i, :])

    # ## solve using casadi
    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp=ocp, solver="ipopt", verbose=False)
    # do some initial iterations to start with a good initial guess
    num_iter_initial = 10
    for _ in range(num_iter_initial):
        casadi_ocp_solver.solve_for_x0(x0_bar = simX_casadi[0,:])
    for i in range(Nsim):
        yref = np.concatenate((X_ref[i, :], U_ref[i, :]))
        for stage in range(ocp.solver_options.N_horizon):
            casadi_ocp_solver.set(stage, "yref", yref)
        casadi_ocp_solver.set(ocp.solver_options.N_horizon, "yref", X_ref[i, :])
        simU_casadi[i, :] = casadi_ocp_solver.solve_for_x0(simX_casadi[i, :])
        simX_casadi[i+1, :] = integrator.simulate(simX_casadi[i, :], simU_casadi[i, :])

    diff_x = np.linalg.norm(simX - simX_casadi)
    diff_u = np.linalg.norm(simU - simU_casadi)
    print(f"Difference in state trajectories: {diff_x}")
    print(f"Difference in control inputs: {diff_u}")
    if diff_x > 1e-4 or diff_u > 1e-4:
        raise Exception("Results with acados and casadi are different!")

    if PLOT:
        Fmax = 80
        acados_u = simU
        acados_x = simX
        casadi_u = simU_casadi
        casadi_x = simX_casadi
        plot_pendulum(np.linspace(0, Tsim, Nsim+1), Fmax, acados_u, acados_x, latexify=False)
        plot_pendulum(np.linspace(0, Tsim, Nsim+1), Fmax, casadi_u, casadi_x, latexify=False)

if __name__ == "__main__":
    main()
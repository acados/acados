#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np
import scipy.linalg
from casadi import vertcat

def setup_ocp_solver(x0, Fmax, N_horizon, Tf):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # set cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    ocp.model.cost_y_expr = vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])

    ocp.constraints.x0 = x0
    ocp.constraints.idxbu = np.array([0])

    # set prediction horizon
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = Tf

    # set ocp options
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.qp_tol = 1e-8

    ocp.code_export_directory = 'c_generated_code_ocp'
    ocp_solver = AcadosOcpSolver(ocp)

    return ocp, ocp_solver


def setup_integrator(ocp: AcadosOcp):
    sim = AcadosSim.from_ocp(ocp)

    sim.solver_options.num_steps = 2 # make integrator more precise than the integrator used within the OCP
    sim.code_export_directory = 'c_generated_code_sim'
    integrator = AcadosSimSolver(sim)
    return integrator


def main():

    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    Fmax = 80

    Tf = .8
    N_horizon = 40

    ocp, ocp_solver = setup_ocp_solver(x0, Fmax, N_horizon, Tf)
    integrator = setup_integrator(ocp)

    nx = ocp_solver.ocp.dims.nx
    nu = ocp_solver.ocp.dims.nu

    Nsim = 100
    simX = np.zeros((Nsim+1, nx))
    simU = np.zeros((Nsim, nu))

    simX[0,:] = x0

    t = np.zeros((Nsim))

    # closed loop
    for i in range(Nsim):

        # solve ocp and get next control input
        simU[i,:] = ocp_solver.solve_for_x0(x0_bar = simX[i, :])

        t[i] = ocp_solver.get_stats('time_tot')

        # simulate system
        simX[i+1, :] = integrator.simulate(x=simX[i, :], u=simU[i,:])

    # evaluate timings
    # scale to milliseconds
    t *= 1000
    print(f'Computation time in ms: min {np.min(t):.3f} median {np.median(t):.3f} max {np.max(t):.3f}')

    # plot results
    model = ocp_solver.ocp.model
    plot_pendulum(np.linspace(0, (Tf/N_horizon)*Nsim, Nsim+1), Fmax, simU, simX, latexify=False, time_label=model.t_label, x_labels=model.x_labels, u_labels=model.u_labels)

    ocp_solver = None


if __name__ == '__main__':
    main()

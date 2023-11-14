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

import sys
sys.path.insert(0, 'common')

from acados_template import AcadosOcp, AcadosSim, AcadosOcpSolver, AcadosSimSolver, latexify_plot
from pendulum_model import export_pendulum_ode_model, export_linearized_pendulum, export_pendulum_ode_model_with_discrete_rk4

from utils import plot_pendulum
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt

X0 = np.array([1.0, np.pi/6, 0.0, 0.0])
FMAX = 80
T_HORIZON = 2.0
N = 40

def create_ocp_solver():
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_linearized_pendulum(0*X0, np.array([0.]))
    # model = export_linearized_pendulum(X0, np.array([1.]))
    # model = export_pendulum_ode_model_with_discrete_rk4(dT=T_HORIZON/N)

    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx

    # set dimensions
    ocp.dims.N = N
    # NOTE: all dimensions but N are now detected automatically in the Python
    #  interface, all other dimensions will be overwritten by the detection.

    # set cost module
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 3*np.diag([1e-1])
    # R = 2*np.diag([1e0])

    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.W_e = Q

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    ocp.constraints.lbu = np.array([-FMAX])
    ocp.constraints.ubu = np.array([+FMAX])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = X0

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    # ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI
    ocp.solver_options.sim_method_num_steps = 1
    ocp.solver_options.tol = 1e-5
    # ocp.solver_options.sim_method_newton_tol = 1e-5

    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.qp_solver_warm_start = 0

    ocp.solver_options.qp_solver_iter_max = 200
    ocp.solver_options.nlp_solver_max_iter = 500
    # ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    ocp.solver_options.nlp_solver_step_length = 1.0
    # set prediction horizon
    ocp.solver_options.tf = T_HORIZON

    acados_ocp_solver = AcadosOcpSolver(ocp)

    return acados_ocp_solver


def create_integrator():
    sim = AcadosSim()
    model = export_pendulum_ode_model()
    sim.model = model
    sim.solver_options.T = T_HORIZON/N
    acados_integrator = AcadosSimSolver(sim)

    return acados_integrator

def main():
    acados_ocp_solver = create_ocp_solver()
    acados_integrator = create_integrator()
    print("Please check the documentation fo eval_param_sens for the requirements on exact solution sensitivities with acados.")

    nx = acados_ocp_solver.acados_ocp.dims.nx
    nu = acados_ocp_solver.acados_ocp.dims.nu

    Nsim = 100
    simX = np.ndarray((Nsim+1, nx))
    simU = np.ndarray((Nsim, nu))
    sens_u = np.ndarray((nu, nx))
    sens_x = np.ndarray((nx, nx))

    xcurrent = X0
    simX[0,:] = xcurrent

    k_lin_feedback = 20 # use lin feedback k_lin_feedback -1 times
    # closed loop
    for i in range(Nsim):
        if i % k_lin_feedback == 0:
            # solve ocp
            acados_ocp_solver.set(0, "lbx", xcurrent)
            acados_ocp_solver.set(0, "ubx", xcurrent)

            status = acados_ocp_solver.solve()

            if status != 0:
                print(xcurrent)
                acados_ocp_solver.print_statistics()
                raise Exception('acados acados_ocp_solver returned status {} in closed loop {}. Exiting.'.format(status, i))

            simU[i,:] = acados_ocp_solver.get(0, "u")

            # calculate solution sensitivities
            u_lin = simU[i,:]
            x_lin = xcurrent

            for index in range(nx):
                acados_ocp_solver.eval_param_sens(index)
                sens_u[:, index] = acados_ocp_solver.get(0, "sens_u")
                sens_x[:, index] = acados_ocp_solver.get(0, "sens_x")
        else:
            # use linear feedback
            # print("using solution sensitivities as feedback")
            simU[i,:] = u_lin + sens_u @ (xcurrent - x_lin)
            # clip
            if simU[i,:] > FMAX:
                simU[i,:] = FMAX
            elif simU[i,:] < -FMAX:
                simU[i,:] = -FMAX
            # for debugging: compute difference of linear feedback wrt real optimal u
            # u_exact = acados_ocp_solver.solve_for_x0(xcurrent)
            # print(f"difference u_lin wrt u_exact {u_exact-simU[i,:]=}")


        # simulate system
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", simU[i,:])

        status = acados_integrator.solve()
        if status != 0:
            raise Exception('acados integrator returned status {}. Exiting.'.format(status))

        # update state
        xcurrent = acados_integrator.get("x")
        simX[i+1,:] = xcurrent

    # plot results
    plot_pendulum(np.linspace(0, T_HORIZON/N*Nsim, Nsim+1), FMAX, simU, simX)


if __name__ == "__main__":
    main()

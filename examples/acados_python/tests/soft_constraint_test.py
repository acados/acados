#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
sys.path.insert(0, '../getting_started/common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver
from export_pendulum_ode_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np
import scipy.linalg

nFormulations = 2
tol = 1E-6

# FORMULATION = 1 # 0 for soft constraints on x - using bounds on x
#                 # 1 for soft constraints on x - using nonlinear constraint h
def run_closed_loop_experiment(FORMULATION):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx
    N = 20

    # set dimensions
    ocp.dims.nx  = nx 
    ocp.dims.ny  = ny 
    ocp.dims.ny_e = ny_e 
    ocp.dims.nu  = model.u.size()[0]
    ocp.dims.ns = nu 
    ocp.dims.N   = N
    ocp.dims.nbu = 1

    if FORMULATION == 0:
        ocp.dims.nh   = 0
        ocp.dims.nsh  = 0
        ocp.dims.nbx  = 1
        ocp.dims.nsbx = 1
    elif FORMULATION == 1:
        ocp.dims.nh   = 1
        ocp.dims.nsh  = 1
        ocp.dims.nbx  = 0
        ocp.dims.nsbx = 0

    # set cost module
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.W_e = Q

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    ocp.cost.zl = 2000*np.ones((1,))
    ocp.cost.Zl = 1*np.ones((1,))
    ocp.cost.zu = 2000*np.ones((1,))
    ocp.cost.Zu = 1*np.ones((1,))

    # set constraints
    Fmax = 80
    vmax = 5

    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0

    # bound on u
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])
    if FORMULATION == 0:
        # soft bound on x
        ocp.constraints.lbx = np.array([-vmax])
        ocp.constraints.ubx = np.array([+vmax])
        ocp.constraints.idxbx = np.array([2]) # v is x[2]
        # indices of slacked constraints within bx
        ocp.constraints.idxsbx = np.array([0])
        # bounds on slack variables
        ocp.constraints.lsbx = np.zeros((ocp.dims.nsbx, ))
        ocp.constraints.usbx = np.zeros((ocp.dims.nsbx, ))

    elif FORMULATION == 1:
        # soft bound on x, using constraint h
        v1 = ocp.model.x[2]
        ocp.model.con_h_expr = v1

        ocp.constraints.lh = np.array([-vmax])
        ocp.constraints.uh = np.array([+vmax])
        # indices of slacked constraints within h
        ocp.constraints.idxsh = np.array([0])
        # bounds on slack variables
        ocp.constraints.lsh = np.zeros((ocp.dims.nh, ))
        ocp.constraints.ush = np.zeros((ocp.dims.nh, ))

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'

    json_filename = 'pendulum_soft_constraints.json'
    acados_ocp_solver = AcadosOcpSolver(ocp, json_file = json_filename)
    acados_integrator = AcadosSimSolver(ocp, json_file = json_filename)

    # closed loop
    Nsim = 20
    simX = np.ndarray((Nsim+1, nx))
    simU = np.ndarray((Nsim, nu))
    xcurrent = x0

    for i in range(Nsim):
        simX[i,:] = xcurrent

        # solve ocp
        acados_ocp_solver.set(0, "lbx", xcurrent)
        acados_ocp_solver.set(0, "ubx", xcurrent)

        status = acados_ocp_solver.solve()
        if status != 0:
            raise Exception('acados acados_ocp_solver returned status {}. Exiting.'.format(status))

        simU[i,:] = acados_ocp_solver.get(0, "u")

        # simulate system
        acados_integrator.set("x", xcurrent)
        acados_integrator.set("u", simU[i,:])

        status = acados_integrator.solve()
        if status != 0:
            raise Exception('acados integrator returned status {}. Exiting.'.format(status))

        # update state
        xcurrent = acados_integrator.get("x")

    simX[Nsim,:] = xcurrent

    # plot results
    plot_pendulum(Tf/N, Fmax, simU, simX, latexify=False)

    # store results
    np.savetxt('test_results/simX_soft_formulation_' + str(FORMULATION), simX)
    np.savetxt('test_results/simU_soft_formulation_' + str(FORMULATION), simU)

    print("soft constraint example: ran formulation", FORMULATION, "successfully.")

if __name__ == "__main__":

    for i in range(nFormulations):
        run_closed_loop_experiment(i)

    simX_ref = np.loadtxt('test_results/simX_soft_formulation_' + str(0))
    simU_ref = np.loadtxt('test_results/simU_soft_formulation_' + str(0))
    for i in range(1, nFormulations):
        simX = np.loadtxt('test_results/simX_soft_formulation_' + str(i))
        simU = np.loadtxt('test_results/simU_soft_formulation_' + str(i))

        error_x = np.linalg.norm(simX_ref - simX)
        error_u = np.linalg.norm(simU_ref - simU)

        error_xu = max([error_x, error_u])

        if error_xu > tol:
            raise Exception("soft constraint example: formulations should return same solution up to" + str(tol))
        else:
            print("soft constraint example: formulation", i, " solution deviates from reference by", error_xu, ".")

    print("soft constraint example: SUCCESS, got same solutions for equivalent formulations up to tolerance %1.e." %(tol))


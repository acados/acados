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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
# from pendulum_model import export_pendulum_ode_model
import numpy as np
# from utils import plot_pendulum
import casadi as cs
import scipy.linalg

def main():
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # parameters
    Tf = 0.15
    N = 3
    nx = 2
    nu = 2
    ny = nx + nu
    ny_e = nx
    x = cs.SX.sym('x', nx)
    u = cs.SX.sym('u', nu)

    A = np.eye(nx)*0.9
    B = np.tril(np.ones((nx, nu)))

    # Initial guess
    x_init = np.outer(np.linspace(1, 0, N+1), [2, 0])
    u_init = np.zeros((N, nu))

    # set model
    model = AcadosModel()
    model.disc_dyn_expr = A @ x + B @ u
    model.x = x
    model.u = u
    # model.p = p
    model.name = "LinearModel"
    ocp.model = model

    # set dimensions
    ocp.dims.N = N

    # set cost
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'
    Q_mat = 0.05*np.diag([1e3, 1e-2])
    R_mat = 0.05*np.diag([1.0, 1.0])

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[ny-1,0] = 1.0
    ocp.cost.Vu = Vu
    ocp.cost.Vx_e = np.eye(nx)
    
    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))
    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    # set constraints
    ocp.constraints.x0 = np.array([2.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' 
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.hessian_approx = 'EXACT' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'DISCRETE'
    # ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP

    # set prediction horizon
    ocp.solver_options.tf = Tf

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_sqp_ocp.json')

    simX_sqp = np.ndarray((N+1, nx))
    simU_sqp = np.ndarray((N, nu))

    status = ocp_solver.solve()
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
    
    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX_sqp[i,:] = ocp_solver.get(i, "x")
        simU_sqp[i,:] = ocp_solver.get(i, "u")
    simX_sqp[N,:] = ocp_solver.get(N, "x")

    print("Solution x: ", simX_sqp)
    print("Solution u: ", simU_sqp)

    ###########################################################################
    # Change options to DDP
    ###########################################################################
    simX_ddp = np.ndarray((N+1, nx))
    simU_ddp = np.ndarray((N, nu))
    ocp.solver_options.nlp_solver_type = 'DDP'

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ddp_ocp.json')

    status = ocp_solver.solve()
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
    
    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX_ddp[i,:] = ocp_solver.get(i, "x")
        simU_ddp[i,:] = ocp_solver.get(i, "u")
    simX_ddp[N,:] = ocp_solver.get(N, "x")

    assert np.allclose(simX_ddp, simX_sqp), "x of ddp and sqp do not coincide"
    assert np.allclose(simU_ddp, simU_sqp), "u of ddp and sqp do not coincide"

    print("Experiment was succesful!")

if __name__ == '__main__':
    main()
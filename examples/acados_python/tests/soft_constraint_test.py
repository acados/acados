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
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np
import scipy.linalg

import itertools

SOFT_CONSTRAINT_TYPES = ['bx', 'h']
TOL = 1E-6
N = 20

QP_SOLVERS = ('PARTIAL_CONDENSING_HPIPM', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM', \
                'FULL_CONDENSING_DAQP'
                )
# no soft constraint support:
# 'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP', \


def run_closed_loop_experiment(soft_constr_type='bx', verbose=False, qp_solver='PARTIAL_CONDENSING_HPIPM'):
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

    # set dimensions
    # NOTE: all dimensions but N ar detected
    ocp.dims.N = N

    # set cost module
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.W_e = Q_mat

    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    Fmax = 80
    vmax = 5

    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0

    # bound on u
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    # soft constraint on x (either via bx or h)
    if soft_constr_type == 'bx':
        # soft bound on x
        ocp.constraints.lbx = np.array([-vmax])
        ocp.constraints.ubx = np.array([+vmax])
        ocp.constraints.idxbx = np.array([2]) # v is x[2]
        # indices of slacked constraints within bx
        ocp.constraints.idxsbx = np.array([0])

    elif soft_constr_type == 'h':
        # soft bound on x, using constraint h
        v1 = ocp.model.x[2]
        ocp.model.con_h_expr = v1

        ocp.constraints.lh = np.array([-vmax])
        ocp.constraints.uh = np.array([+vmax])
        # indices of slacked constraints within h
        ocp.constraints.idxsh = np.array([0])
    else:
        raise Exception(f"soft_constr_type must be 'bx', or 'h', got {soft_constr_type}.")

    ocp.cost.zl = 2000*np.ones((1,))
    ocp.cost.Zl = 1*np.ones((1,))
    ocp.cost.zu = 2000*np.ones((1,))
    ocp.cost.Zu = 1*np.ones((1,))

    # set options
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.tol = 1e-1 * TOL
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    ocp.solver_options.qp_solver_warm_start = 0
    ocp.solver_options.qp_solver_iter_max = 200

    json_filename = 'pendulum_soft_constraints.json'
    acados_ocp_solver = AcadosOcpSolver(ocp, json_file = json_filename)
    acados_integrator = AcadosSimSolver(ocp, json_file = json_filename)

    # closed loop
    Nsim = 20
    simX = np.ndarray((Nsim+1, nx))
    simU = np.ndarray((Nsim, nu))
    xcurrent = x0
    qp_iter = np.zeros(Nsim)
    sqp_iter = np.zeros(Nsim)

    for i in range(Nsim):
        simX[i,:] = xcurrent

        # solve ocp
        acados_ocp_solver.set(0, "lbx", xcurrent)
        acados_ocp_solver.set(0, "ubx", xcurrent)

        status = acados_ocp_solver.solve()
        if verbose:
            acados_ocp_solver.print_statistics()

        if status != 0:
            acados_ocp_solver.print_statistics()
            raise Exception('acados acados_ocp_solver returned status {}. Exiting.'.format(status))

        qp_iter[i] = np.sum(acados_ocp_solver.get_stats('qp_iter'))
        sqp_iter[i] = acados_ocp_solver.get_stats('sqp_iter')

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

    # get slack values at stage 1
    sl = acados_ocp_solver.get(1, "sl")
    su = acados_ocp_solver.get(1, "su")
    print("sl", sl, "su", su)

    # plot results
    plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)

    # store results
    np.savetxt(f'test_results/simX_soft_formulation_{soft_constr_type}_{qp_solver}', simX)
    np.savetxt(f'test_results/simU_soft_formulation_{soft_constr_type}_{qp_solver}', simU)
    np.savetxt(f'test_results/sqp_iter_soft_formulation_{soft_constr_type}_{qp_solver}', sqp_iter)

    print(f"\nsoft constraint example: ran formulation {soft_constr_type} with {qp_solver} successfully.\n")
    print(f"took {sqp_iter} SQP iterations and {qp_iter} QP iterations.")

    del acados_ocp_solver

def main():
    # import pdb; pdb.set_trace()
    for (soft_constr_type, qp_solver) in itertools.product(SOFT_CONSTRAINT_TYPES, QP_SOLVERS):
        run_closed_loop_experiment(soft_constr_type=soft_constr_type, qp_solver=qp_solver)

    # load reference solution
    soft_constr_type = 'bx'
    qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    simX_ref = np.loadtxt(f'test_results/simX_soft_formulation_{soft_constr_type}_{qp_solver}')
    simU_ref = np.loadtxt(f'test_results/simU_soft_formulation_{soft_constr_type}_{qp_solver}')
    sqp_iter_ref = np.loadtxt(f'test_results/sqp_iter_soft_formulation_{soft_constr_type}_{qp_solver}')

    # compare
    for (soft_constr_type, qp_solver) in itertools.product(SOFT_CONSTRAINT_TYPES, QP_SOLVERS):
        simX = np.loadtxt(f'test_results/simX_soft_formulation_{soft_constr_type}_{qp_solver}')
        simU = np.loadtxt(f'test_results/simU_soft_formulation_{soft_constr_type}_{qp_solver}')
        sqp_iter = np.loadtxt(f'test_results/sqp_iter_soft_formulation_{soft_constr_type}_{qp_solver}')

        error_x = np.linalg.norm(simX_ref - simX)
        error_u = np.linalg.norm(simU_ref - simU)

        error_xu = max([error_x, error_u])

        print(f"soft constraint example: formulation {soft_constr_type} with {qp_solver} deviates from reference by {error_xu}")

        if error_xu > TOL:
            raise Exception(f"soft constraint example: formulations should return same solution up to {TOL:.2e}, got error_x {error_x}, error_u {error_u} for {soft_constr_type}, {qp_solver}")

        if any(sqp_iter != sqp_iter_ref):
            raise Exception(f"all formulations should take the same number of SQP iterations.")

    print(f"soft constraint example: SUCCESS, got same solutions for equivalent formulations up to tolerance {TOL:.2e}")



if __name__ == "__main__":
    main()
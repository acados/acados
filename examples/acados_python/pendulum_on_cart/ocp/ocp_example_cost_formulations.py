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
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosMultiphaseOcp
from pendulum_model import export_pendulum_ode_model, export_augmented_pendulum_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
import casadi as ca
from casadi.tools import entry, struct_symSX

COST_VERSIONS = ['LS', 'EXTERNAL', 'EXTERNAL_Z', 'NLS', 'NLS_TO_EXTERNAL', 'NLS_Z', 'LS_Z', 'CONL', 'CONL_Z', 'AUTO']
HESSIAN_APPROXIMATION = 'GAUSS_NEWTON' # 'GAUSS_NEWTON
N = 20
T_HORIZON = 1.0

NX = 4
NU = 1
FMAX = 80

def formulate_ocp(cost_version: str, constraint_version="bu") -> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    if cost_version in ['EXTERNAL_Z','NLS_Z', 'LS_Z', 'CONL_Z']:
        model = export_augmented_pendulum_model()
        nz = model.z.rows()
    else:
        model = export_pendulum_ode_model()

    model.name += cost_version

    # set model
    ocp.model = model

    ny = NX + NU
    ny_e = NX

    ocp.solver_options.N_horizon = N

    # set cost
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    x = ocp.model.x
    u = ocp.model.u

    cost_W = scipy.linalg.block_diag(Q_mat, R_mat)

    if cost_version in ['LS', 'NLS', 'NLS_TO_EXTERNAL', 'NLS_Z', 'LS_Z', 'CONL', 'CONL_Z']:
        ocp.cost.yref = np.zeros((ny, ))
        ocp.cost.yref_e = np.zeros((ny_e, ))
    if cost_version in ['LS', 'NLS', 'NLS_TO_EXTERNAL', 'NLS_Z', 'LS_Z']:
        ocp.cost.W_e = Q_mat
        ocp.cost.W = cost_W

    if cost_version in ['CONL', 'CONL_Z', 'EXTERNAL', 'EXTERNAL_Z', 'AUTO']:
        cost_W = ca.sparsify(ca.DM(cost_W))
        Q_mat = ca.sparsify(ca.DM(Q_mat))

    if cost_version == 'LS':
        ocp.cost.cost_type = 'LINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'

        ocp.cost.Vx = np.zeros((ny, NX))
        ocp.cost.Vx[:NX,:NX] = np.eye(NX)

        Vu = np.zeros((ny, NU))
        Vu[4,0] = 1.0
        ocp.cost.Vu = Vu

        ocp.cost.Vx_e = np.eye(NX)
        ocp.solver_options.fixed_hess = 1

    elif cost_version == 'LS_Z':
        ocp.cost.cost_type = 'LINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'

        ocp.cost.Vx = np.zeros((ny, NX))
        ocp.cost.Vx[:NX,:NX] = np.eye(NX)
        ocp.cost.Vx[0, 0] = 0.0

        ocp.cost.Vz = np.zeros((ny, nz))
        ocp.cost.Vz[0, 0] = 1.0

        Vu = np.zeros((ny, NU))
        Vu[4,0] = 1.0
        ocp.cost.Vu = Vu
        ocp.cost.Vx_e = np.eye(NX)

    elif cost_version == 'NLS':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'

        ocp.model.cost_y_expr = ca.vertcat(x, u)
        ocp.model.cost_y_expr_e = x

    elif cost_version == 'NLS_Z':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'

        y_expr_z = ca.vertcat(model.z[0], model.x[1:], u)
        print(f"{y_expr_z}")
        ocp.model.cost_y_expr = y_expr_z
        ocp.model.cost_y_expr_e = x

    elif cost_version == 'CONL':
        ocp.cost.cost_type = 'CONVEX_OVER_NONLINEAR'
        ocp.cost.cost_type_e = 'CONVEX_OVER_NONLINEAR'

        ocp.model.cost_y_expr = ca.vertcat(x, u)
        ocp.model.cost_y_expr_e = x

        r = ca.SX.sym('r', ny)
        r_e = ca.SX.sym('r_e', ny_e)
        ocp.model.cost_r_in_psi_expr = r
        ocp.model.cost_r_in_psi_expr_e = r_e

        ocp.model.cost_psi_expr = 0.5 * (r.T @ cost_W @ r)
        ocp.model.cost_psi_expr_e = 0.5 * (r_e.T @ Q_mat @ r_e)

        # with custom hessian
        # ocp.model.cost_conl_custom_outer_hess = cost_W

    elif cost_version == 'CONL_Z':
        ocp.cost.cost_type = 'CONVEX_OVER_NONLINEAR'
        ocp.cost.cost_type_e = 'CONVEX_OVER_NONLINEAR'

        r = ca.SX.sym('r', ny)
        r_e = ca.SX.sym('r_e', ny_e)

        ocp.model.cost_psi_expr = 0.5 * (r.T @ cost_W @ r)
        ocp.model.cost_psi_expr_e = 0.5 * (r_e.T @ Q_mat @ r_e)

        ocp.model.cost_r_in_psi_expr = r
        ocp.model.cost_r_in_psi_expr_e = r_e

        y_expr_z = ca.vertcat(model.z[0], model.x[1:], u)
        print(f"{y_expr_z}")
        ocp.model.cost_y_expr = y_expr_z
        ocp.model.cost_y_expr_e = x

    elif cost_version == 'EXTERNAL':
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'

        ocp.model.cost_expr_ext_cost = .5*ca.vertcat(x, u).T @ cost_W @ ca.vertcat(x, u)
        ocp.model.cost_expr_ext_cost_e = .5*x.T @ Q_mat @ x

    elif cost_version == 'EXTERNAL_Z':
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'

        y_expr_z = ca.vertcat(model.z[0], model.x[1:], u)

        ocp.model.cost_expr_ext_cost = .5*y_expr_z.T @ cost_W @ y_expr_z
        ocp.model.cost_expr_ext_cost_e = .5*x.T @ Q_mat @ x

    elif cost_version == 'AUTO':
        ocp.cost.cost_type = 'AUTO'
        ocp.cost.cost_type_e = 'AUTO'
        ocp.model.cost_expr_ext_cost = .5*ca.vertcat(x, u).T @ cost_W @ ca.vertcat(x, u)
        ocp.model.cost_expr_ext_cost_e = .5*x.T @ Q_mat @ x

    elif cost_version == 'NLS_TO_EXTERNAL':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'

        ocp.model.cost_y_expr = ca.vertcat(x, u)
        ocp.model.cost_y_expr_e = x
        ocp.translate_cost_to_external_cost(cost_hessian='GAUSS_NEWTON')
    elif cost_version == 'NLS_TO_EXTERNAL_P_GLOBAL':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'

        p_global = struct_symSX([
            entry('W', shape=(ny, ny)),
            entry('yref', shape=(ny, )),
            entry('W_e', shape=(ny_e, ny_e)),
            entry('yref_e', shape=(ny_e, ))
        ])
        ocp.model.p_global = p_global.cat

        ocp.cost.W = p_global['W']
        ocp.cost.yref = p_global['yref']
        ocp.cost.W_e = p_global['W_e']
        ocp.cost.yref_e = p_global['yref_e']

        ocp.model.cost_y_expr = ca.vertcat(x, u)
        ocp.model.cost_y_expr_e = x

        ocp.translate_cost_to_external_cost(cost_hessian='GAUSS_NEWTON')

        p_global_values = p_global(0)
        p_global_values['W'] = cost_W
        p_global_values['yref'] = np.zeros((ny, ))
        p_global_values['W_e'] = Q_mat
        p_global_values['yref_e'] = np.zeros((ny_e, ))

        ocp.p_global_values = p_global_values.cat.full().flatten()
    else:
        raise Exception('Unknown cost_version.')

    # set constraints
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
    if constraint_version == "bu":
        ocp.constraints.lbu = np.array([-FMAX])
        ocp.constraints.ubu = np.array([+FMAX])
        ocp.constraints.idxbu = np.array([0])
    elif constraint_version == "h":
        ocp.constraints.lh = np.array([-FMAX])
        ocp.constraints.uh = np.array([+FMAX])
        ocp.model.con_h_expr = ocp.model.u[0]  # h = u[0]
        ocp.constraints.lh_0 = np.array([-FMAX])
        ocp.constraints.uh_0 = np.array([+FMAX])
        ocp.model.con_h_expr_0 = ocp.model.u[0]  # h = u[0]

    return ocp


def main(cost_version: str, formulation_type='ocp', integrator_type='IRK', reformulate_to_external=False, plot=False):

    if cost_version == 'EXTERNAL':
        ext_cost_use_num_hess = True
    else:
        ext_cost_use_num_hess = False

    if formulation_type == 'mocp':
        ocp = AcadosMultiphaseOcp(N_list=[1, N-1])

        phase_0 = formulate_ocp(cost_version)
        phase_1 = formulate_ocp(cost_version)
        ocp.set_phase(phase_0, 0)
        ocp.set_phase(phase_1, 1)
        ocp.solver_options = ocp.solver_options
        ocp.name = 'mocp_' + phase_0.model.name
        if isinstance(integrator_type, list):
            ocp.mocp_opts.integrator_type = integrator_type
        ocp.solver_options.sim_method_num_steps = np.array([1] + (N-1)*[5])
    else:
        ocp = formulate_ocp(cost_version)

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES, FULL_CONDENSING_DAQP, FULL_CONDENSING_HPIPM
    ocp.solver_options.hessian_approx = HESSIAN_APPROXIMATION

    if not isinstance(integrator_type, list):
        ocp.solver_options.integrator_type = integrator_type

    ocp.solver_options.qp_solver_cond_N = 5

    # set prediction horizon
    ocp.solver_options.tf = T_HORIZON
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI

    if cost_version in ['EXTERNAL', 'EXTERNAL_Z']:
        ocp.solver_options.ext_cost_num_hess = ext_cost_use_num_hess
    # from casadi import jacobian
    # ux = ca.vertcat(ocp.model.u, ocp.model.x)
    # jacobian(jacobian(ocp.model.cost_expr_ext_cost, ux), ux)
    # ca.SX(@1=0.04, @2=4000,
    # [[@1, 00, 00, 00, 00],
    #  [00, @2, 00, 00, 00],
    #  [00, 00, @2, 00, 00],
    #  [00, 00, 00, @1, 00],
    #  [00, 00, 00, 00, @1]])

    if reformulate_to_external:
        ocp.solver_options.fixed_hess = 0
        if cost_version in ['LS', 'NLS', 'NLS_Z', 'LS_Z', 'CONL', 'CONL_Z']:
            p = ca.SX.sym('yref', ocp.cost.yref.shape[0])
            p_values = ocp.cost.yref
            yref = p
            yref_e = p[:ocp.cost.yref_e.shape[0]]
        else:
            p = p_values = yref= yref_e = None
        ocp.translate_cost_to_external_cost(p=p, p_values=p_values, yref=yref, yref_e=yref_e)

    # create solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    # NOTE: hessian is wrt [u,x]
    if ext_cost_use_num_hess and cost_version in  ['EXTERNAL', 'EXTERNAL_Z']:
        for i in range(N):
            ocp_solver.cost_set(i, "ext_cost_num_hess", np.diag([0.02, 2000, 2000, 0.02, 0.02, ]))
        ocp_solver.cost_set(N, "ext_cost_num_hess", np.diag([2000, 2000, 0.02, 0.02, ]))

    # test set cost scaling
    if formulation_type == 'ocp':
        for i in range(N):
            ocp_solver.cost_set(i, "scaling", ocp.solver_options.time_steps[i])

    simX = np.zeros((N+1, NX))
    simU = np.zeros((N, NU))

    status = ocp_solver.solve()

    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N,:] = ocp_solver.get(N, "x")

    cost_val = ocp_solver.get_cost()

    # compare with reference solution
    cost_val_ref = 2242.166372615175
    rel_diff_cost = abs(cost_val - cost_val_ref) / cost_val_ref
    print(f"cost value is: {cost_val}, difference with reference: {rel_diff_cost:.2e}")
    if rel_diff_cost > 1e-6:
        raise Exception(f"Cost value is not correct: rel_diff_cost = {rel_diff_cost:.2e} > {1e-6}.")

    # test getter
    if formulation_type == 'mocp':
        cost = ocp.cost[0]
        cost_e = ocp.cost[-1]
    else:
        cost = cost_e = ocp.cost

    if cost.cost_type in ['LINEAR_LS', 'NONLINEAR_LS', 'CONVEX_OVER_NONLINEAR']:
        yref_ = ocp_solver.cost_get(1, 'yref')
        assert np.allclose(yref_, cost.yref)

    if cost.cost_type in ['LINEAR_LS', 'NONLINEAR_LS']:
        W_ = ocp_solver.cost_get(1, 'W')
        assert np.allclose(W_, cost.W)

    if cost_e.cost_type_e in ['LINEAR_LS', 'NONLINEAR_LS', 'CONVEX_OVER_NONLINEAR']:
        yref_e_ = ocp_solver.cost_get(ocp.solver_options.N_horizon, 'yref')
        assert np.allclose(yref_e_, cost_e.yref_e)

    if cost.cost_type in ['LINEAR_LS', 'NONLINEAR_LS']:
        W_e_ = ocp_solver.cost_get(ocp.solver_options.N_horizon, 'W')
        assert np.allclose(W_e_, cost_e.W_e)

    # plot results
    if plot:
        plot_pendulum(np.linspace(0, T_HORIZON, N+1), FMAX, simU, simX, latexify=False)

    ocp_solver.store_iterate(filename='solution.json', overwrite=True)
    ocp_solver.load_iterate(filename='solution.json')

if __name__ == "__main__":
    main(cost_version='LS', formulation_type='mocp', integrator_type=['IRK', 'ERK'], plot=False)
    for cost_version in COST_VERSIONS:
        for formulation_type in ['ocp', 'mocp']:
            print(f"cost version: {cost_version}, formulation type: {formulation_type}")
            main(cost_version=cost_version, formulation_type=formulation_type, plot=False)

    for cost_version in ["NLS_TO_EXTERNAL_P_GLOBAL"]:
        print(f"cost version: {cost_version} reformulated as EXTERNAL cost")
        main(cost_version=cost_version, formulation_type='ocp', plot=False, reformulate_to_external=True)

# timings
# time_tot = 1e8
# time_lin = 1e8

# for k in range(1000):

#     status = ocp_solver.solve()
#     time_tot = min(time_tot, ocp_solver.get_stats("time_tot")[0])
#     time_lin = min(time_lin, ocp_solver.get_stats("time_lin")[0])

# print("CPU time = ", time_tot * 1e3, "ms")
# print("CPU time linearization = ", time_lin * 1e3, "ms")
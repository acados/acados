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


import casadi as ca
from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY
import numpy as np

def get_hs099_definition():
    # The optimal objective is (if given in):
    f_opt = -0.831079892e+9
    x_opt = ca.DM([0.542468, 0.529021, 0.508449, 0.480269, 0.451236, 0.409183, 0.352788])
    x = ca.MX.sym('x', 7)
    x0 = np.zeros((7, 1))
    lbx = -ACADOS_INFTY*np.ones((7, 1))
    ubx = ACADOS_INFTY*np.ones((7, 1))
    a = ca.DM.zeros(8, 1)
    t = ca.DM.zeros(8, 1)
    r_tmp = ca.MX.zeros(8, 1)
    q_tmp = ca.MX.zeros(8, 1)
    s_tmp = ca.MX.zeros(8, 1)
    g = ca.MX.zeros(2, 1)
    lbg = -ACADOS_INFTY*np.ones((2, 1))
    ubg = ACADOS_INFTY*np.ones((2, 1))

    a[0] = 0
    a[1] = 50
    a[2] = 50
    a[3] = 75
    a[4] = 75
    a[5] = 75
    a[6] = 100
    a[7] = 100
    b = 32
    t[0] = 0
    t[1] = 25
    t[2] = 50
    t[3] = 100
    t[4] = 150
    t[5] = 200
    t[6] = 290
    t[7] = 380

    x0[0:7] = 0.5
    lbx[0:7] = 0
    ubx[0:7] = 1.58

    r_tmp[0] = 0
    for i in range(1, 8):
        r_tmp[i] = a[i]*(t[i]-t[i-1]) * ca.cos(x[i-1]) + r_tmp[i-1]

    s_tmp[0] = 0
    for i in range(1, 8):
        s_tmp[i] = (t[i]-t[i-1])*(a[i]*ca.sin(x[i-1]) - b) + s_tmp[i-1]
    q_tmp[0] = 0
    for i in range(1, 8):
        q_tmp[i] = 0.5*(t[i]-t[i-1])**2*(a[i]*ca.sin(x[i-1]) - b) + (t[i]-t[i-1])*s_tmp[i-1] + q_tmp[i-1]

    lbg[0] = 1.0e+5
    ubg[0] = 1.0e+5
    g[0] = q_tmp[7]
    lbg[1] = 1.0e+3
    ubg[1] = 1.0e+3
    g[1] = s_tmp[7]

    f = -r_tmp[7]**2

    return (x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)


def create_hs099_ocp(x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()
    ocp.code_export_directory = 'c_generated_code_hs099'

    # set model
    model = ocp.model
    model.x = x
    model.name = 'hs099'
    ocp.cost.cost_type_e = 'EXTERNAL'
    model.cost_expr_ext_cost_e = f
    model.con_h_expr_e = g

    ocp.constraints.lh_e = lbg
    ocp.constraints.uh_e = ubg

    ocp.constraints.idxbx_e = np.arange(7)
    ocp.constraints.lbx_e = lbx
    ocp.constraints.ubx_e = ubx

    ocp.solver_options.N_horizon = 0

    return ocp


def set_ocp_options(ocp: AcadosOcp, use_qp_scaling: bool = True, qp_tol_strategy: str = 'NAIVE'):
    opts = ocp.solver_options
    opts.hessian_approx = 'EXACT'
    opts.nlp_solver_max_iter = 20
    opts.qp_solver_iter_max = 100

    opts.nlp_solver_ext_qp_res = 1
    opts.nlp_solver_type = 'SQP'

    opts.globalization = 'FIXED_STEP'

    # Scaling
    if use_qp_scaling:
        opts.qpscaling_scale_objective = 'OBJECTIVE_GERSHGORIN'
        opts.qpscaling_scale_constraints = 'INF_NORM'
        opts.qpscaling_lb_norm_inf_grad_obj = 1e-4
        opts.qpscaling_ub_max_abs_eig = 1e5
    if qp_tol_strategy == 'NAIVE':
        pass
    elif qp_tol_strategy == 'SUFFICIENTLY_SMALL':
        opts.qp_tol = 1e-9
    elif qp_tol_strategy == 'ADAPTIVE':
        opts.nlp_qp_tol_strategy = 'ADAPTIVE_QPSCALING'
    return

def test_solver(use_qp_scaling: bool = True, qp_tol_strategy: str = 'NAIVE'):

    x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0 = get_hs099_definition()

    ocp = create_hs099_ocp(x_opt, f_opt, x, f, g, lbg, ubg, lbx, ubx, x0)
    set_ocp_options(ocp, use_qp_scaling, qp_tol_strategy)
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    # initialize
    ocp_solver.set(0, 'x', x0)

    # solve
    status = ocp_solver.solve()

    # evaluate
    ocp_solver.print_statistics()

    x_sol = ocp_solver.get(0, 'x')

    diff_x = np.linalg.norm(x_sol - x_opt)
    if diff_x > 1e-6:
        raise ValueError(f"Solution x does not match optimal x: {diff_x}")

    cost = ocp_solver.get_cost()
    if not np.isclose(cost, f_opt, atol=1e-6):
        raise ValueError(f"Cost does not match optimal cost: {cost} vs {f_opt}")

    if status != 0:
        raise ValueError(f"Solver failed with status: {status}")


def main():
    try:
        test_solver(use_qp_scaling=True, qp_tol_strategy='NAIVE')
    except Exception as e:
        if 'status' in str(e):
            print(f"Test failed with status error as expected: {e}")
        else:
            raise ValueError(f"Test with QP scaling and NAIVE strategy should fail")
    else:
        raise ValueError("Test with QP scaling should not fail!")

    test_solver(use_qp_scaling=True, qp_tol_strategy='SUFFICIENTLY_SMALL')
    test_solver(use_qp_scaling=True, qp_tol_strategy='ADAPTIVE')
    test_solver(use_qp_scaling=False)

if __name__ == "__main__":
    main()

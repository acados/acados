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

from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY, AcadosOcpFlattenedIterate
import numpy as np
import casadi as ca


def create_solver(solver_name: str, nlp_solver_type: str = 'SQP_WITH_FEASIBLE_QP',
                  allow_switching_modes: bool = True,
                  use_qp_scaling: bool = False,
                  soft_h: bool = True):

    ocp = AcadosOcp()

    nx = 2
    # set model
    ocp.model.name = f"dense_nlp_{solver_name}"
    ocp.model.x = ca.SX.sym('x', nx, 1)

    ny = nx

    # discretization
    N = 0

    ocp.cost.W_e = 2*np.diag([1e3, 1e3])

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx_e = np.eye((nx))
    ocp.cost.yref_e = np.ones((ny, ))

    # set constraints
    xmax = 2.0
    ocp.constraints.lbx_e = -xmax * np.ones((nx,))
    ocp.constraints.ubx_e = +xmax * np.ones((nx,))
    ocp.constraints.idxbx_e = np.arange(nx)
    # soften the bounds
    # ocp.constraints.idxsbx_e = np.array([0, 1])
    # ocp.cost.zl_e = np.array([1.0, 1.0])
    # ocp.cost.zu_e = np.array([1e0, 1e0])
    # ocp.cost.Zl_e = np.array([1.0, 1.0])
    # ocp.cost.Zu_e = np.array([1e0, 1e0])

    # define soft nonlinear constraint
    scale_h = 1.0
    radius = 1.0
    ocp.model.con_h_expr_e = scale_h * (ocp.model.x[0]**2 + ocp.model.x[1]**2)
    ocp.constraints.lh_e = -1000 * np.ones((1,))
    ocp.constraints.lh_e = -ACADOS_INFTY * np.ones((1,))
    ocp.constraints.uh_e = scale_h * radius**2 * np.ones((1,))

    # soften
    if soft_h:
        ocp.constraints.idxsh_e = np.array([0])
        ocp.cost.zl_e = np.concatenate((ocp.cost.zl_e, np.array([1.0])))
        ocp.cost.zu_e = np.concatenate((ocp.cost.zu_e, np.array([1e0])))
        ocp.cost.Zl_e = np.concatenate((ocp.cost.Zl_e, np.array([1.0])))
        ocp.cost.Zu_e = np.concatenate((ocp.cost.Zu_e, np.array([1e0])))

    # set options
    solver_options = ocp.solver_options
    solver_options.N_horizon = N

    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    qp_tol = 5e-9
    solver_options.qp_tol = qp_tol
    solver_options.qp_solver_ric_alg = 1
    solver_options.qp_solver_mu0 = 1e4
    solver_options.qp_solver_iter_max = 400
    solver_options.hessian_approx = 'GAUSS_NEWTON'
    solver_options.nlp_solver_type = nlp_solver_type
    solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    solver_options.globalization_full_step_dual = True
    solver_options.print_level = 1
    solver_options.nlp_solver_max_iter = 6
    solver_options.use_constraint_hessian_in_feas_qp = False
    solver_options.nlp_solver_ext_qp_res = 0

    if not allow_switching_modes:
        solver_options.search_direction_mode = 'BYRD_OMOJOKUN'
        solver_options.allow_direction_mode_switch_to_nominal = False

    if use_qp_scaling:
        ocp.solver_options.qpscaling_scale_constraints = "INF_NORM"
        ocp.solver_options.qpscaling_scale_objective = "OBJECTIVE_GERSHGORIN"

    # create ocp solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    return ocp, ocp_solver

def check_qp_scaling(ocp_solver: AcadosOcpSolver):
    qpscaling_status = ocp_solver.get_stats("qpscaling_status")
    if qpscaling_status == 0:
        print("QP scaling reported no issues.")
    else:
        print(f"QP scaling reported issues with status {qpscaling_status}.")

    if ocp_solver.acados_ocp.solver_options.qpscaling_scale_constraints == "NO_CONSTRAINT_SCALING":
        try:
            constraint_scaling = ocp_solver.get_qp_scaling_constraints(0)
        except Exception as e:
            print(f"constraint scaling not done as expected.")
            return

    for i in range(ocp_solver.N+1):
        constraint_scaling = ocp_solver.get_qp_scaling_constraints(i)
        print(f"Constraint scaling at stage {i}: {constraint_scaling}")
        # if not np.all(constraint_scaling != 1.0):
        #     raise ValueError(f"Constraint scaling should have non-unit to actually test the functionality")
    objective_scaling = ocp_solver.get_qp_scaling_objective()
    print(f"Objective scaling: {objective_scaling}")


def call_solver(ocp_solver: AcadosOcpSolver) -> AcadosOcpFlattenedIterate:
    # solve
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    sqp_iter = ocp_solver.get_stats('sqp_iter')
    if status != 0:
        # raise RuntimeError(f"acados returned status {status} after {sqp_iter} SQP iterations.")
        print(f'acados returned status {status}.')

    print(f"cost function value = {ocp_solver.get_cost()} after {sqp_iter} SQP iterations")
    sol = ocp_solver.store_iterate_to_flat_obj()
    return sol

def check_solutions(sol_1: AcadosOcpFlattenedIterate, sol_2: AcadosOcpFlattenedIterate, soft_h: bool):
    # check solutions
    for field in ["x", "u", "sl", "su", "lam", "pi"]:
        v1 = getattr(sol_1, field)
        v2 = getattr(sol_2, field)
        if not np.allclose(v1, v2, atol=1e-6):
            print(f"Field {field} differs: max diff = {np.max(np.abs(v1 - v2))}")
            print(f"got difference {v1 - v2}")
        else:
            print(f"Solutions match in field {field}.")
            pass
    print(f"{sol_1}")

    if soft_h:
        if np.any(sol_1.su > 1e-1):
            print("checked with active soft constraints.")
        else:
            raise ValueError("Soft constraints should be active, but are not.")

    if sol_1.allclose(sol_2):
        print("Both solvers have the same solution.")
    else:
        raise ValueError("Solutions of solvers differ!")

def check_residual_solutions(stat1: np.ndarray, stat2: np.ndarray):
    n_rows1 = len(stat1[0])
    n_rows2 = len(stat2[0])

    assert n_rows1 == n_rows2, f"Both solvers should take the same number of iterations!, got {n_rows1} for solver 1, and {n_rows2} for solver 2"

    for jj in range(n_rows1):
        # res_stat
        assert np.allclose(stat1[1][jj], stat2[1][jj]), f"res_stat differs in iter {jj}"
        # res_eq
        assert np.allclose(stat1[2][jj], stat2[2][jj]), f"res_eq differs in iter {jj}"
        # res_ineq
        assert np.allclose(stat1[3][jj], stat2[3][jj]), f"res_ineq differs in iter {jj}"
        # res_comp
        assert np.allclose(stat1[4][jj], stat2[4][jj]), f"res_comp differs in iter {jj}"

def test_qp_scaling(soft_h: bool = True):
    nlp_solver_type = "SQP"
    # nlp_solver_type = "SQP_WITH_FEASIBLE_QP"

    # test without QP scaling
    print("Reference ...")
    _, ocp_solver_2 = create_solver("2", nlp_solver_type=nlp_solver_type, allow_switching_modes=True, use_qp_scaling=False, soft_h=soft_h)
    sol_2 = call_solver(ocp_solver_2)
    stats2 = ocp_solver_2.get_stats("statistics")

    check_qp_scaling(ocp_solver_2)
    print(f"Reference solution: {sol_2}")

    # test QP scaling
    print("Testing QP scaling with SQP solver...")
    _, ocp_solver_1 = create_solver("1", nlp_solver_type=nlp_solver_type, allow_switching_modes=True, use_qp_scaling=True, soft_h=soft_h)
    sol_1 = call_solver(ocp_solver_1)
    check_qp_scaling(ocp_solver_1)
    stats1 = ocp_solver_1.get_stats("statistics")

    check_residual_solutions(stats1, stats2)
    check_solutions(sol_1, sol_2, soft_h)

def test_sanity_check(soft_h: bool = True, use_qp_scaling: bool = True):
    print("Sanity Check SQP and SQP_WITH_FEASIBLE_QP solver...")

    # test without QP scaling
    print("Solving with SQP")
    _, ocp_solver_2 = create_solver("2", nlp_solver_type="SQP", allow_switching_modes=True, use_qp_scaling=use_qp_scaling, soft_h=soft_h)
    sol_2 = call_solver(ocp_solver_2)
    check_qp_scaling(ocp_solver_2)
    stats2 = ocp_solver_2.get_stats("statistics")
    print(f"Reference solution: {sol_2}")

    # test QP scaling
    print("Solving with SQP_WITH_FEASIBLE_QP")
    _, ocp_solver_1 = create_solver("1", nlp_solver_type="SQP_WITH_FEASIBLE_QP", allow_switching_modes=True, use_qp_scaling=use_qp_scaling, soft_h=soft_h)
    sol_1 = call_solver(ocp_solver_1)
    stats1 = ocp_solver_1.get_stats("statistics")
    check_qp_scaling(ocp_solver_1)

    check_residual_solutions(stats1, stats2)
    check_solutions(sol_1, sol_2, soft_h)
    print("\n")


if __name__ == '__main__':
    # Sanity Checks
    test_sanity_check(soft_h=False, use_qp_scaling=False) # overall solver sanity check
    test_sanity_check(soft_h=True, use_qp_scaling=False) # overall solver sanity check
    test_sanity_check(soft_h=False, use_qp_scaling=True)
    test_sanity_check(soft_h=True, use_qp_scaling=True)

    test_qp_scaling(soft_h=False)
    test_qp_scaling(soft_h=True)



#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

# Binary constraints with PARTIAL_CONDENSING_DAQP: branch and bound vs. brute force over all
# fixed binary assignments (each a QP solved by FULL_CONDENSING_DAQP).

import itertools
import numpy as np
import casadi as ca
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, ACADOS_INFTY

N = 6
# binary: input s in {0, 1} (u), the same as a general constraint (g), or velocity v in {-0.1, 0.1} (x)
VALUES = {'u': (0.0, 1.0), 'g': (0.0, 1.0), 'x': (-0.1, 0.1)}


def create_ocp(variant: str, qp_solver: str) -> AcadosOcp:
    # double integrator x = [p, v], u = [a, s]; the binary variable is listed first in idxbu/idxbx
    model = AcadosModel()
    model.name = f'binary_{variant}_{qp_solver[:4].lower()}'
    p, v, a, s = ca.SX.sym('p'), ca.SX.sym('v'), ca.SX.sym('a'), ca.SX.sym('s')
    model.x = ca.vertcat(p, v)
    model.u = ca.vertcat(a, s)
    model.disc_dyn_expr = model.x + 0.1 * ca.vertcat(v, a + s)

    ocp = AcadosOcp()
    ocp.model = model
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = 0.1 * N
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'

    ocp.cost.cost_type = ocp.cost.cost_type_e = 'LINEAR_LS'
    ocp.cost.Vx = np.vstack([np.eye(2), np.zeros((2, 2))])
    ocp.cost.Vu = np.vstack([np.zeros((2, 2)), np.eye(2)])
    ocp.cost.W = np.diag([10, 1, 0.1, 0.1])
    ocp.cost.yref = np.array([0, 0, 0, 0.3])
    ocp.cost.Vx_e = np.eye(2)
    ocp.cost.W_e = np.diag([10, 1])
    ocp.cost.yref_e = np.zeros(2)

    c = ocp.constraints
    c.x0 = np.zeros(2)
    binary = qp_solver == 'PARTIAL_CONDENSING_DAQP'
    if variant == 'u':
        c.idxbu, c.lbu, c.ubu = np.array([1, 0]), np.array([0, -1]), np.array([1, 1])
        if binary:
            c.idxbinu = np.array([1])
    elif variant == 'g':
        c.idxbu, c.lbu, c.ubu = np.array([0]), np.array([-1]), np.array([1])
        c.C, c.D, c.lg, c.ug = np.zeros((1, 2)), np.array([[0.0, 1.0]]), np.array([0.0]), np.array([1.0])
        if binary:
            c.idxbing = np.array([0])
    else:
        c.idxbu, c.lbu, c.ubu = np.array([0, 1]), np.array([-1, -1]), np.array([1, 1])
        c.idxbx, c.lbx, c.ubx = np.array([1, 0]), np.array([-0.1, -5]), np.array([0.1, 5])
        if binary:
            c.idxbinx = np.array([1])
    ocp.code_gen_options.code_export_directory = f'c_generated_code_{model.name}'
    ocp.code_gen_options.json_file = f'{model.name}.json'
    return ocp


def lam_sum(solver: AcadosOcpSolver):
    return sum(np.sum(np.abs(solver.get(t, 'lam'))) for t in range(N + 1))


def solve(solver: AcadosOcpSolver, variant, x0, lower, upper):
    """Solve with the given range of the binary variable at each stage; returns its values"""
    stages = range(1, N) if variant == 'x' else range(N)
    for t, lo, up in zip(stages, lower, upper):
        if variant == 'u':
            solver.constraints_set(t, 'lbu', np.array([lo, -1]))
            solver.constraints_set(t, 'ubu', np.array([up, 1]))
        elif variant == 'g':
            solver.constraints_set(t, 'lg', np.array([lo]))
            solver.constraints_set(t, 'ug', np.array([up]))
        else:
            solver.constraints_set(t, 'lbx', np.array([lo, -5]))
            solver.constraints_set(t, 'ubx', np.array([up, 5]))
    solver.reset()
    solver.set(0, 'lbx', x0)
    solver.set(0, 'ubx', x0)
    status = solver.solve()
    z = [solver.get(t, 'x')[1] if variant == 'x' else solver.get(t, 'u')[1] for t in stages]
    return status, solver.get_cost(), np.array(z)


def test_binary_constraints(variant):
    solver = AcadosOcpSolver(create_ocp(variant, 'PARTIAL_CONDENSING_DAQP'), verbose=False)
    reference = AcadosOcpSolver(create_ocp(variant, 'FULL_CONDENSING_DAQP'), verbose=False)
    lo, up = VALUES[variant]
    nbin = N - 1 if variant == 'x' else N
    is_binary = lambda z: np.minimum(np.abs(z - lo), np.abs(z - up)) < 1e-8
    H = 2  # binary_horizon: binaries at the stages t < H (x: stages 1, ..., H-1)
    nbin_H = H - 1 if variant == 'x' else H
    fractional = False

    for x0 in [np.array([1.0, 0.0]), np.array([-0.3, 0.2]), np.array([0.05, 0.0]), np.array([0.0, -0.1])]:
        status, cost, z = solve(solver, variant, x0, [lo] * nbin, [up] * nbin)
        assert status == 0, f'{variant}, x0 = {x0}: status {status}'
        assert np.all(is_binary(z)), f'{variant}, x0 = {x0}: not binary: {z}'
        assert lam_sum(solver) == 0, f'{variant}, x0 = {x0}: multipliers after branch and bound'

        results = [solve(reference, variant, x0, combo, combo) for combo in itertools.product((lo, up), repeat=nbin)]
        best = min(c for st, c, _ in results if st == 0)
        assert abs(cost - best) <= 1e-6 * max(1, best), f'{variant}, x0 = {x0}: cost {cost}, brute force {best}'

        # fixed binaries: same cost
        status, cost_fixed, _ = solve(solver, variant, x0, z, z)
        assert status == 0 and abs(cost_fixed - cost) <= 1e-6 * max(1, cost), \
            f'{variant}, x0 = {x0}: fixed: {cost_fixed} != {cost}'
        set_bounds = ([lo] * nbin, [up] * nbin)

        # relaxed binaries: a lower bound, with multipliers; then binary again
        solver.options_set('qp_relax_binaries', 1)
        status, cost_relaxed, z_relaxed = solve(solver, variant, x0, *set_bounds)
        assert status == 0 and cost_relaxed <= cost + 1e-8, f'{variant}, x0 = {x0}: relaxed {cost_relaxed} > {cost}'
        assert lam_sum(solver) > 0, f'{variant}, x0 = {x0}: no multipliers with relaxed binaries'
        fractional |= not np.all(is_binary(z_relaxed))
        solver.options_set('qp_relax_binaries', 0)
        status, cost_again, z_again = solve(solver, variant, x0, *set_bounds)
        assert status == 0 and abs(cost_again - cost) <= 1e-6 * max(1, cost) and np.all(is_binary(z_again)), \
            f'{variant}, x0 = {x0}: binary again: {cost_again} != {cost}'

        # binary horizon: binaries at the first stages only, the cost between the relaxed and the binary one
        solver.options_set('qp_binary_horizon', H)
        status, cost_H, z_H = solve(solver, variant, x0, *set_bounds)
        solver.options_set('qp_binary_horizon', 0)
        assert status == 0 and np.all(is_binary(z_H[:nbin_H])), f'{variant}, x0 = {x0}: binary horizon: {z_H}'
        assert cost_relaxed - 1e-8 <= cost_H <= cost + 1e-8, \
            f'{variant}, x0 = {x0}: binary horizon: {cost_H} not in [{cost_relaxed}, {cost}]'
        print(f'{variant}, x0 = {x0}: cost {cost:.6f}, brute force {best:.6f}, binary horizon {H} {cost_H:.6f}, '
              f'relaxed {cost_relaxed:.6f}, binaries {np.round(z, 6) + 0}')
    assert fractional, f'{variant}: the relaxed binaries are never fractional'

    # the options of binary constraints with another QP solver: an error (not an exit of the QP solver)
    for field in ['qp_relax_binaries', 'qp_binary_horizon']:
        try:
            reference.options_set(field, 1)
        except ValueError:
            continue
        raise AssertionError(f'expected a ValueError for {field} with FULL_CONDENSING_DAQP')


def test_binary_errors():
    def expect_error(modify):
        ocp = create_ocp('u', 'PARTIAL_CONDENSING_DAQP')
        modify(ocp)
        try:
            ocp.make_consistent()
        except ValueError:
            return
        raise AssertionError('expected a ValueError')

    expect_error(lambda ocp: setattr(ocp.constraints, 'idxbinu', np.array([2])))  # no bound
    expect_error(lambda ocp: setattr(ocp.constraints, 'idxbing', np.array([0])))  # no general constraint
    expect_error(lambda ocp: setattr(ocp.constraints, 'lbu', np.array([-ACADOS_INFTY, -1])))
    expect_error(lambda ocp: setattr(ocp.solver_options, 'qp_solver', 'PARTIAL_CONDENSING_HPIPM'))
    expect_error(lambda ocp: setattr(ocp.solver_options, 'qp_solver_cond_N', N // 2))


if __name__ == '__main__':
    test_binary_errors()
    for variant in ['u', 'g', 'x']:
        test_binary_constraints(variant)

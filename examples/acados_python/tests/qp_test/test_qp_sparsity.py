#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcpQp, AcadosOcpQpSolver, AcadosOcpQpOptions, AcadosOcpIterate
import numpy as np
import pytest

N_HORIZON = 20
SPARSITY_FIELDS = ['A', 'B', 'Q', 'R', 'S', 'C', 'D']
# solvers that use the sparsity pattern
SOLVERS = ['PARTIAL_CONDENSING_OSQP']
N_TIMING_SOLVES = 5
# longer chains show larger speedups, but take minutes: 30 masses and more are left out here
N_MASS_VALUES = [3, 10, 20]
# the three ways a QP is solved by each solver
VARIANTS = ['no sparsity', 'default pattern', 'pattern of the data']


def create_chain_qp(n_mass: int, soft: bool) -> AcadosOcpQp:
    """
    Free chain of n_mass masses connected by springs, with damping and one force input per mass.
    The states are the positions and the velocities. The dynamics are discretized with one explicit Euler step,
    such that A and B keep the sparsity of the chain.
    The inputs are bounded and the distance of neighboring masses is limited by general constraints,
    which are soft if requested.
    """
    N = N_HORIZON
    nx, nu, ng = 2 * n_mass, n_mass, n_mass - 1
    dt, stiffness, damping = 0.1, 2.0, 0.1
    u_max, max_distance = 0.5, 0.15

    spring = -2.0 * np.eye(n_mass) + np.eye(n_mass, k=1) + np.eye(n_mass, k=-1)
    spring[0, 0] = -1.0
    spring[-1, -1] = -1.0
    A = np.block([[np.eye(n_mass), dt * np.eye(n_mass)],
                  [dt * stiffness * spring, (1.0 - dt * damping) * np.eye(n_mass)]])
    B = np.vstack([np.zeros((n_mass, n_mass)), dt * np.eye(n_mass)])
    Q = np.diag(np.concatenate([10.0 * np.ones(n_mass), 0.1 * np.ones(n_mass)]))
    R = 0.1 * np.eye(nu)
    # distance of neighboring masses
    C = np.hstack([np.eye(ng, n_mass, k=1) - np.eye(ng, n_mass), np.zeros((ng, n_mass))])
    # the masses start 0.1 apart, moving in alternating directions
    x0 = np.concatenate([0.1 * np.arange(n_mass), 0.1 * (-1.0) ** np.arange(n_mass)])

    qp = AcadosOcpQp(N)
    for stage in range(N + 1):
        nu_stage = nu if stage < N else 0
        ng_stage = ng if stage > 0 else 0
        nbx = nx if stage == 0 else 0
        nb = nu_stage + nbx

        # every stage gets its own copy of the matrices
        qp.set('Q', stage, Q.copy())
        qp.set('q', stage, np.zeros(nx))
        qp.set('R', stage, R[:nu_stage, :nu_stage].copy())
        qp.set('S', stage, np.zeros((nu_stage, nx)))
        qp.set('r', stage, np.zeros(nu_stage))
        if stage < N:
            qp.set('A', stage, A.copy())
            qp.set('B', stage, B.copy())
            qp.set('b', stage, np.zeros(nx))

        # input bounds, and the initial state as equality constrained state bounds
        qp.set('idxb', stage, np.arange(nb))
        qp.set('lbu', stage, -u_max * np.ones(nu_stage))
        qp.set('ubu', stage, u_max * np.ones(nu_stage))
        qp.set('lbx', stage, x0[:nbx])
        qp.set('ubx', stage, x0[:nbx])
        if stage == 0:
            qp.set('idxe', stage, np.arange(nu_stage, nb))

        qp.set('C', stage, C[:ng_stage, :].copy())
        qp.set('D', stage, np.zeros((ng_stage, nu_stage)))
        qp.set('lg', stage, -max_distance * np.ones(ng_stage))
        qp.set('ug', stage, max_distance * np.ones(ng_stage))

        ns = ng_stage if soft else 0
        idxs_rev = -np.ones(nb + ng_stage)
        if soft:
            idxs_rev[nb:] = np.arange(ns)
        qp.set('idxs_rev', stage, idxs_rev)
        for field in ['Zl', 'Zu']:
            qp.set(field, stage, 100.0 * np.ones(ns))
        for field in ['zl', 'zu']:
            qp.set(field, stage, 10.0 * np.ones(ns))
        qp.set('lls', stage, np.zeros(ns))
        qp.set('lus', stage, np.zeros(ns))

        for field, n in [('lbu', nu_stage), ('ubu', nu_stage), ('lbx', nbx), ('ubx', nbx),
                         ('lg', ng_stage), ('ug', ng_stage), ('lls', ns), ('lus', ns)]:
            qp.set(f'{field}_mask', stage, np.ones(n))
    return qp


def create_solver(qp: AcadosOcpQp, qp_solver: str, cond_N: int, with_sparsity: bool = False, tol: float = 1e-8) -> AcadosOcpQpSolver:
    opts = AcadosOcpQpOptions()
    opts.qp_solver = qp_solver
    opts.cond_N = cond_N
    opts.iter_max = 20000
    opts.tol_stat = opts.tol_eq = opts.tol_ineq = opts.tol_comp = tol
    # every solve starts from zero, such that repeated solves of the same QP can be timed.
    # NOTE: OSQP keeps the step size it adapted in earlier solves, the iteration count settles after a few solves.
    opts.warm_start = 0
    opts.with_sparsity = with_sparsity
    return AcadosOcpQpSolver(qp, opts=opts)


def get_sparsity_of_data(qp: AcadosOcpQp, stage: int) -> dict:
    # the pattern of a matrix: 1 where the matrix is nonzero
    fields = [field for field in SPARSITY_FIELDS if stage < qp.N or field not in ['A', 'B']]
    return {field: (getattr(qp, field)[stage] != 0.0).astype(float) for field in fields}


def set_sparsity_of_data(solver: AcadosOcpQpSolver, qp: AcadosOcpQp):
    for stage in range(qp.N + 1):
        for field, pattern in get_sparsity_of_data(qp, stage).items():
            solver.set(stage, f'{field}_sparsity', pattern)


def count_entries(qp: AcadosOcpQp):
    n_entries, n_nonzero = 0, 0
    for stage in range(qp.N + 1):
        for pattern in get_sparsity_of_data(qp, stage).values():
            n_entries += pattern.size
            n_nonzero += int(np.sum(pattern))
    return n_entries, n_nonzero


def solve(solver: AcadosOcpQpSolver, case: str) -> AcadosOcpIterate:
    status = solver.solve()
    assert status == 0, f"{case}: {solver.qp_solver_name} returned status {status}"
    return solver.get_iterate()


def max_difference(sol: AcadosOcpIterate, ref: AcadosOcpIterate) -> float:
    return max(np.max(np.abs(a - b)) for field in ['x', 'u', 'sl', 'su']
               for a, b in zip(getattr(sol, field), getattr(ref, field)) if a.size > 0)


def median_solve_time(solver: AcadosOcpQpSolver) -> float:
    # total time of a solve: condensing, update of the solver data, solver call and expansion
    times = []
    for _ in range(N_TIMING_SOLVES):
        solver.solve()
        times.append(solver.get_stats('time_tot'))
    return float(np.median(times))


def run_same_solution_case(qp_solver: str, n_mass: int, soft: bool, cond_N: int):
    """
    Solves the QP
    a) without sparsity patterns,
    b) with the default sparsity patterns, which keep all entries,
    c) with the sparsity patterns of the data,
    and checks that the solutions are the same. Prints the solve times, next to the one of HPIPM.
    """
    case = f"n_mass = {n_mass}, soft = {soft}, cond_N = {cond_N}"
    qp = create_chain_qp(n_mass, soft)

    ref = solve(create_solver(qp, 'PARTIAL_CONDENSING_HPIPM', cond_N, tol=1e-10), case)
    scale = max(np.max(np.abs(u)) for u in ref.u if u.size > 0)
    assert scale > 0.1, f"{case}: the test problem should have a nonzero solution"

    solvers = {
        'no sparsity': create_solver(qp, qp_solver, cond_N),
        'default pattern': create_solver(qp, qp_solver, cond_N, with_sparsity=True),
        'pattern of the data': create_solver(qp, qp_solver, cond_N, with_sparsity=True),
    }
    set_sparsity_of_data(solvers['pattern of the data'], qp)

    # same solution in all three variants, and the same as the reference
    solutions = {name: solve(solver, f"{case}, {name}") for name, solver in solvers.items()}
    diff_ref = {}
    for name, sol in solutions.items():
        diff = max_difference(sol, solutions['no sparsity'])
        assert diff < 1e-6, f"{case}: solution with {name} differs by {diff:.2e} from the one without sparsity for {qp_solver}"
        diff_ref[name] = max_difference(sol, ref)
        assert diff_ref[name] < 1e-4, f"{case}: solution with {name} differs by {diff_ref[name]:.2e} from the reference for {qp_solver}"

    times = {name: median_solve_time(solver) for name, solver in solvers.items()}
    time_hpipm = median_solve_time(create_solver(qp, 'PARTIAL_CONDENSING_HPIPM', cond_N))

    # the solutions are still the same after solving again with the same matrix structure
    for name, solver in solvers.items():
        diff = max_difference(solver.get_iterate(), solutions['no sparsity'])
        assert diff < 1e-6, f"{case}: solution with {name} differs by {diff:.2e} after solving again for {qp_solver}"

    n_entries, n_nonzero = count_entries(qp)
    speedup = times['no sparsity'] / times['pattern of the data']
    print(f"{n_mass:>6} {str(soft):>5} {cond_N:>6} | {n_entries:>7} {n_nonzero:>7} | {1e3*time_hpipm:>10.3f} | "
          + " | ".join(f"{1e3*times[name]:>{len(name) + 5}.3f}" for name in VARIANTS)
          + f" | {speedup:>7.2f} | {max(diff_ref.values()):>23.1e}")


def run_entry_outside_pattern_case(qp_solver: str):
    """An entry that is nonzero in the data, but zero in the pattern, is treated as zero."""
    n_mass, cond_N, stage = 4, N_HORIZON, 2
    qp = create_chain_qp(n_mass, soft=False)
    ref = solve(create_solver(qp, 'PARTIAL_CONDENSING_HPIPM', cond_N, tol=1e-10), "reference")

    # the spring between the first two masses, as seen by the first mass
    i, j = n_mass, 1
    assert qp.A[stage][i, j] != 0.0
    A_zeroed = qp.A[stage].copy()
    A_zeroed[i, j] = 0.0
    qp_zeroed = create_chain_qp(n_mass, soft=False)
    qp_zeroed.set('A', stage, A_zeroed)
    ref_zeroed = solve(create_solver(qp_zeroed, 'PARTIAL_CONDENSING_HPIPM', cond_N, tol=1e-10), "reference with zeroed entry")
    assert max_difference(ref_zeroed, ref) > 1e-3, "the zeroed entry should change the solution"

    solver = create_solver(qp, qp_solver, cond_N, with_sparsity=True)
    set_sparsity_of_data(solver, qp)
    pattern = get_sparsity_of_data(qp, stage)['A']
    pattern[i, j] = 0.0
    solver.set(stage, 'A_sparsity', pattern)
    sol = solve(solver, "entry outside the pattern")

    diff = max_difference(sol, ref_zeroed)
    assert diff < 1e-4, f"entry outside the pattern: solution differs by {diff:.2e} from the QP with that entry set to zero for {qp_solver}"
    print(f"entry outside the pattern: {qp_solver} solves the QP with that entry set to zero.")


@pytest.mark.parametrize("qp_solver", SOLVERS)
@pytest.mark.parametrize("n_mass", N_MASS_VALUES)
@pytest.mark.parametrize("soft", [False, True])
@pytest.mark.parametrize("cond_N", [N_HORIZON, N_HORIZON // 4])
def test_qp_sparsity_preserves_solution(qp_solver, n_mass, soft, cond_N):
    run_same_solution_case(qp_solver, n_mass, soft, cond_N)


@pytest.mark.parametrize("qp_solver", SOLVERS)
def test_qp_sparsity_ignores_entry_outside_pattern(qp_solver):
    run_entry_outside_pattern_case(qp_solver)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))

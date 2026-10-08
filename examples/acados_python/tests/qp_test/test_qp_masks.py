#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcpQp, AcadosOcpQpSolver, AcadosOcpQpOptions, AcadosOcpIterate
import numpy as np

QP_FILE = 'last_qp_one_sided_test.json'

# tolerance on x, u, sl, su w.r.t. the reference solution
SOLVER_TOLS = {
    'PARTIAL_CONDENSING_HPIPM': 1e-6,
    'PARTIAL_CONDENSING_OSQP': 1e-4,
    'PARTIAL_CONDENSING_CLARABEL': 1e-4,
}
SOLVERS_WITH_SLACK_BOUND_MASKS = ['PARTIAL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_OSQP', 'PARTIAL_CONDENSING_CLARABEL']


def create_solver(qp: AcadosOcpQp, qp_solver: str, warm_start: int = 0) -> AcadosOcpQpSolver:
    opts = AcadosOcpQpOptions()
    opts.qp_solver = qp_solver
    opts.iter_max = 4000
    opts.warm_start = warm_start
    return AcadosOcpQpSolver(qp, opts=opts)


def solve_reference(qp: AcadosOcpQp) -> AcadosOcpIterate:
    solver = create_solver(qp, 'PARTIAL_CONDENSING_HPIPM')
    status = solver.solve()
    assert status == 0, f"reference solver returned status {status}"
    return solver.get_iterate()


def get_mask(qp: AcadosOcpQp, stage: int) -> np.ndarray:
    # same ordering as lam
    return np.concatenate([qp.lbu_mask[stage], qp.lbx_mask[stage], qp.lg_mask[stage],
                           qp.ubu_mask[stage], qp.ubx_mask[stage], qp.ug_mask[stage],
                           qp.lls_mask[stage], qp.lus_mask[stage]])


def check_solution(solver: AcadosOcpQpSolver, qp: AcadosOcpQp, ref: AcadosOcpIterate, case: str):
    qp_solver = solver.qp_solver_name
    assert solver.status == 0, f"{case}: {qp_solver} returned status {solver.status}"
    sol = solver.get_iterate()
    for stage in range(qp.N + 1):
        for field in ['x', 'u', 'sl', 'su']:
            np.testing.assert_allclose(getattr(sol, field)[stage], getattr(ref, field)[stage], atol=SOLVER_TOLS[qp_solver],
                                       err_msg=f"{case}: {field} mismatch at stage {stage} for {qp_solver}")
        lam = sol.lam[stage]
        mask = get_mask(qp, stage)
        assert np.all(lam[mask == 0.0] == 0.0), f"{case}: nonzero multiplier of masked constraint at stage {stage} for {qp_solver}: lam = {lam}"


def poison_masked_lbx(qp: AcadosOcpQp, ref: AcadosOcpIterate):
    # write values into masked lower state bounds that cut off the reference solution if enforced;
    # only in the second half of the horizon, such that the QP stays feasible if they are enforced
    for stage in range(qp.N // 2, qp.N + 1):
        nu = qp.dims.nu[stage]
        nbu = qp.dims.nbu[stage]
        lbx = qp.lbx[stage].copy()
        for j in np.where(qp.lbx_mask[stage] == 0.0)[0]:
            lbx[j] = ref.x[stage][qp.idxb[stage][nbu + j] - nu] + 0.5
        qp.set('lbx', stage, lbx)


def create_soft_qp(mask_lls: bool) -> AcadosOcpQp:
    # add one slack to the masked lower state bound at stages 1, ..., N-1
    qp = AcadosOcpQp.from_json(QP_FILE)
    for stage in range(1, qp.N):
        qp.set('idxs_rev', stage, np.array([-1, 0]))
        for field in ['Zl', 'Zu', 'zl', 'zu']:
            qp.set(field, stage, np.ones(1))
        qp.set('lls', stage, np.zeros(1))
        qp.set('lus', stage, np.zeros(1))
        qp.set('lls_mask', stage, np.zeros(1) if mask_lls else np.ones(1))
        qp.set('lus_mask', stage, np.ones(1))
    return qp


def run_poisoned_case(qp_solver: str, warm_start: int, ref: AcadosOcpIterate):
    qp = AcadosOcpQp.from_json(QP_FILE)
    poison_masked_lbx(qp, ref)
    solver = create_solver(qp, qp_solver, warm_start)
    solver.solve()
    check_solution(solver, qp, ref, 'poisoned')


def run_flip_case(qp_solver: str, warm_start: int, ref: AcadosOcpIterate):
    # mask an active upper input bound on an existing solver, poison it such that lbu > ubu, restore it
    qp = AcadosOcpQp.from_json(QP_FILE)
    stage = 1
    i_ubu = qp.dims.nb[stage] + qp.dims.ng[stage]
    assert ref.lam[stage][i_ubu] > 1e-3, "ubu has to be active in the reference solution"

    qp_masked = AcadosOcpQp.from_json(QP_FILE)
    qp_masked.set('ubu_mask', stage, np.zeros(1))
    ref_masked = solve_reference(qp_masked)

    solver = create_solver(qp, qp_solver, warm_start)
    solver.solve()
    check_solution(solver, qp, ref, 'flip, original')

    solver.set(stage, 'ubu_mask', np.zeros(1))
    solver.set(stage, 'ubu', qp.lbu[stage] - 1.0)
    solver.solve()
    check_solution(solver, qp_masked, ref_masked, 'flip, masked')

    solver.set(stage, 'ubu_mask', np.ones(1))
    solver.set(stage, 'ubu', qp.ubu[stage])
    solver.solve()
    check_solution(solver, qp, ref, 'flip, restored')


def run_soft_case(qp_solver: str, warm_start: int):
    mask_lls = qp_solver in SOLVERS_WITH_SLACK_BOUND_MASKS
    ref = solve_reference(create_soft_qp(mask_lls))

    qp = create_soft_qp(mask_lls)
    poison_masked_lbx(qp, ref)
    if mask_lls:
        for stage in range(1, qp.N):
            qp.set('lls', stage, np.ones(1))
    solver = create_solver(qp, qp_solver, warm_start)
    solver.solve()
    check_solution(solver, qp, ref, 'soft')


if __name__ == "__main__":
    ref = solve_reference(AcadosOcpQp.from_json(QP_FILE))
    for qp_solver in SOLVER_TOLS:
        for warm_start in [0, 1]:
            print(f"qp_solver = {qp_solver}, warm_start = {warm_start}")
            run_poisoned_case(qp_solver, warm_start, ref)
            run_flip_case(qp_solver, warm_start, ref)
            run_soft_case(qp_solver, warm_start)
    print("test_qp_masks: all checks passed.")

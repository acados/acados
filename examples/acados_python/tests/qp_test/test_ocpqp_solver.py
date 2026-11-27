from acados_template import AcadosOcpQp, AcadosOcpQpSolver, AcadosOcpIterate, AcadosOcpQpOptions
from pathlib import Path

import numpy as np
import pytest

TEST_DIR = Path(__file__).resolve().parent
SQP_QP_SOL_PAIRS = [
    ('last_qp_nonuniform_pendulum.json', 'sqp_sol_nonuniform_pendulum.json'),
    ('last_qp_one_sided_test.json', 'sqp_sol_one_sided_test.json'),
]
QP_SOLVERS = [
    'PARTIAL_CONDENSING_HPIPM',
    'FULL_CONDENSING_HPIPM',
    'FULL_CONDENSING_DAQP',
    'PARTIAL_CONDENSING_OSQP',
    'PARTIAL_CONDENSING_CLARABEL',
    'FULL_CONDENSING_QPOASES',
]

# NOTE: create test cases in SQP example using:
    # id = 'xxx'
    # ocp_solver.dump_last_qp_to_json(filename=f'last_qp_{id}.json', overwrite=True)
    # ocp_solver.store_iterate(filename=f'sqp_sol_{id}.json', overwrite=True)


@pytest.mark.parametrize("qp_solver", QP_SOLVERS)
@pytest.mark.parametrize(("qp_json_file", "sqp_sol_file"), SQP_QP_SOL_PAIRS)
def test_ocpqp_solver_matches_sqp_solution(qp_solver, qp_json_file, sqp_sol_file):
    qp = AcadosOcpQp.from_json(TEST_DIR / qp_json_file)
    opts = AcadosOcpQpOptions()
    opts.iter_max = 500
    opts.qp_solver = qp_solver
    solver = AcadosOcpQpSolver(qp, opts=opts)

    status = solver.solve()
    assert status == 0, f"{qp_solver} returned status {status}"
    sol = solver.get_iterate()

    sqp_sol = AcadosOcpIterate.from_json(TEST_DIR / sqp_sol_file)
    tol = 1e-5

    for stage in range(qp.N + 1):
        np.testing.assert_allclose(
            sol.lam[stage],
            sqp_sol.lam[stage],
            atol=tol,
            err_msg=f"lam mismatch at stage {stage} for {qp_solver}",
        )

    for stage in range(qp.N):
        np.testing.assert_allclose(
            sol.pi[stage],
            sqp_sol.pi[stage],
            atol=tol,
            err_msg=f"pi mismatch at stage {stage} for {qp_solver}",
        )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))

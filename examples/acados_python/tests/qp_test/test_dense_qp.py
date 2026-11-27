from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import numpy as np
import pytest


@pytest.mark.parametrize(
    "solver_name",
    [
        "IPOPT",
        "PARTIAL_CONDENSING_HPIPM",
        "FULL_CONDENSING_HPIPM",
    ],
)
def test_dense_qp(solver_name):
    nv = 2

    qp = AcadosOcpQp(N=0)
    qp.set('Q', 0, 1000*np.eye(nv))
    qp.set('q', 0, np.ones((nv,)))
    qp.set('idxs_rev', 0, -1 * np.ones((nv,)))
    qp.set('lbx', 0, -np.ones((nv,)))
    qp.set('ubx', 0, np.ones((nv,)))
    qp.set('lbx_mask', 0, -np.ones((nv,)))
    qp.set('ubx_mask', 0, np.ones((nv,)))
    qp.set('idxb', 0, np.arange(nv))

    if solver_name == 'IPOPT':
        solver = AcadosCasadiOcpQpSolver(qp)
    else:
        opts = AcadosOcpQpOptions()
        opts.qp_solver = solver_name
        opts.print_level = 0
        solver = AcadosOcpQpSolver(qp, opts)

    status = solver.solve()
    assert status == 0, f"{solver_name} returned status {status}"
    sol = solver.get_iterate()
    expected_solution = -0.001 * np.ones(nv)
    np.testing.assert_allclose(sol.x[0], expected_solution, atol=1e-6)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
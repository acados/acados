from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import numpy as np


def main(solver_name: str = 'HPIPM'):
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

    if solver_name in ['IPOPT']:
        solver = AcadosCasadiOcpQpSolver(qp)
    elif solver_name == 'FULL_CONDENSING_HPIPM':
        opts = AcadosOcpQpOptions()
        opts.qp_solver = 'FULL_CONDENSING_HPIPM'
        opts.print_level = 0
        solver = AcadosOcpQpSolver(qp, opts)
    elif solver_name == 'PARTIAL_CONDENSING_HPIPM':
        opts = AcadosOcpQpOptions()
        opts.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
        opts.print_level = 0
        solver = AcadosOcpQpSolver(qp, opts)
    else:
        raise NotImplementedError(f"solver {solver_name} not available.")

    solver.solve()
    sol = solver.get_iterate()
    print(sol)
    if 'HPIPM' in solver_name:
        solver.print_statistics()

if __name__ == "__main__":
    main(solver_name='IPOPT')
    main(solver_name='PARTIAL_CONDENSING_HPIPM')
    main(solver_name='FULL_CONDENSING_HPIPM')
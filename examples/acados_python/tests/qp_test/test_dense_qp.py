from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import numpy as np


def main(Acasadi:bool = False):
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

    if Acasadi:
        solver = AcadosCasadiOcpQpSolver(qp)
    else:
        opts = AcadosOcpQpOptions()
        opts.print_level = 1
        solver = AcadosOcpQpSolver(qp, opts)

    solver.solve()
    sol = solver.get_iterate()
    print(sol)

if __name__ == "__main__":
    main(Acasadi=True)
    main(Acasadi=False)
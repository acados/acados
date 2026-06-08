from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import numpy as np


def main():
    qp = AcadosOcpQp(N=0)
    nv = 2

    # dummy
    qp.set('R', 0, np.eye(0))
    qp.set('r', 0, np.ones((0,)))
    qp.set('S', 0, np.zeros((0, nv)))
    # TODO: can we make those optional?
    qp.set('lbu', 0, -np.ones((0,)))
    qp.set('ubu', 0, np.ones((0,)))
    qp.set('lls', 0, -np.ones((0,)))
    qp.set('lus', 0, np.ones((0,)))
    qp.set('lbu_mask', 0, -np.ones((0,)))
    qp.set('ubu_mask', 0, np.ones((0,)))
    qp.set('lls_mask', 0, -np.ones((0,)))
    qp.set('lus_mask', 0, np.ones((0,)))
    qp.set('idxe', 0, np.ones(0))
    
    # slack cost
    qp.set('zl', 0, -np.ones((0,)))
    qp.set('zu', 0, np.ones((0,)))
    qp.set('Zl', 0, -np.ones((0,)))
    qp.set('Zu', 0, np.ones((0,)))
    # g empty
    qp.set('C', 0, np.eye(0,0))
    qp.set('D', 0, np.eye(0,0))
    qp.set('lg', 0, -np.ones((0,)))
    qp.set('ug', 0, np.ones((0,)))
    qp.set('lg_mask', 0, -np.ones((0,)))
    qp.set('ug_mask', 0, np.ones((0,)))
    # actual problem
    qp.set('idxs_rev', 0, -1 * np.ones((nv,)))
    qp.set('Q', 0, 1000*np.eye(nv))
    qp.set('q', 0, np.ones((nv,)))
    qp.set('lbx', 0, -np.ones((nv,)))
    qp.set('ubx', 0, np.ones((nv,)))
    qp.set('lbx_mask', 0, -np.ones((nv,)))
    qp.set('ubx_mask', 0, np.ones((nv,)))
    qp.set('idxb', 0, np.arange(nv))

    # TODO: make consistent should not be needed and called inside solver creation.
    qp.make_consistent()
    opts = AcadosOcpQpOptions()
    # TODO: print_level does not work yet at creation time, only via opts_set.
    opts.print_level = 10
    solver = AcadosOcpQpSolver(qp)
    solver.opts_set('print_level', 2)

    solver.solve()
    sol = solver.get_iterate()
    print(sol)
    breakpoint()

if __name__ == "__main__":
    main()
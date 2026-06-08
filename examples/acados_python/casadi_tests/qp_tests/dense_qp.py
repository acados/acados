from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import numpy as np


def main():
    nv = 2
    dense_qp_dict = {
        'N': 0,
        'Q':[1000*np.eye(nv)],
        'q':[np.ones((nv,))],
        'idxs_rev':[-1 * np.ones((nv,))],
        'lbx': [-np.ones((nv,))],
        'ubx': [np.ones((nv,))],
        'lbx_mask': [-np.ones((nv,))],
        'ubx_mask': [np.ones((nv,))],
        'idxb': [np.arange(nv)],
    }

    qp = AcadosOcpQp.from_qp_dict(dense_qp_dict)
    opts = AcadosOcpQpOptions()
    opts.print_level = 1
    solver = AcadosOcpQpSolver(qp, opts)

    solver.solve()
    sol = solver.get_iterate()
    print(sol)

if __name__ == "__main__":
    main()
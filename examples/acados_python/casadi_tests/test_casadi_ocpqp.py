import sys
sys.path.insert(0, '../tests/qp_test')

from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import time
import numpy as np

sqp_qp_sol_pairs = [
        ('last_qp_nonuniform_pendulum.json', 'sqp_sol_nonuniform_pendulum.json'),
                    ]

# NOTE: create test cases in SQP example using:
    # id = 'xxx'
    # ocp_solver.dump_last_qp_to_json(filename=f'last_qp_{id}.json', overwrite=True)
    # ocp_solver.store_iterate(filename=f'sqp_sol_{id}.json', overwrite=True)

for qp_solver in ['PARTIAL_CONDENSING_HPIPM']:
    n_iter_list = []
    timings_list = []
    for qp_json_file, sqp_sol_file in sqp_qp_sol_pairs:
        qp = AcadosOcpQp.from_json(qp_json_file)
        opts = AcadosOcpQpOptions()
        opts.iter_max = 500
        opts.qp_solver = qp_solver

        acados_solver = AcadosOcpQpSolver(qp, opts=opts)
        status = acados_solver.solve()
        iterate = acados_solver.get_iterate()
        acados_u = np.array([acados_solver.get(i, "u") for i in range(qp.N)])
        acados_x = np.array([acados_solver.get(i, "x") for i in range(qp.N+1)])
        acados_lam = np.concatenate([acados_solver.get(i, "lam") for i in range(qp.N+1)])

        casadi_solver = AcadosCasadiOcpQpSolver(qp)
        casadi_solver.set_iterate(iterate)
        status = casadi_solver.solve()
        casadi_u = np.array([casadi_solver.get(i, "u") for i in range(qp.N)])
        casadi_x = np.array([casadi_solver.get(i, "x") for i in range(qp.N+1)])
        casadi_lam = np.concatenate([casadi_solver.get(i, "lam") for i in range(qp.N+1)])

        # evaluate difference
        assert np.allclose(casadi_u, acados_u, atol=1e-5), f"u mismatch for {qp_json_file}"
        assert np.allclose(casadi_x, acados_x, atol=1e-5), f"x mismatch for {qp_json_file}"
        assert np.allclose(casadi_lam, acados_lam, atol=1e-5), f"lam mismatch for {qp_json_file}"
        print(f"QP solution matches for {qp_json_file}.")
from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import numpy as np


qp_json_files = ['examples/acados_python/casadi_tests/qp_tests/prob_0.json',
                'examples/acados_python/casadi_tests/qp_tests/pendulum_qp.json',]
def main():
    for qp_json_file in qp_json_files:
        qp = AcadosOcpQp.from_json(qp_json_file)
        opts = AcadosOcpQpOptions()
        opts.iter_max = 500
        opts.qp_solver = 'PARTIAL_CONDENSING_HPIPM'

        acados_solver = AcadosOcpQpSolver(qp, opts=opts)
        status = acados_solver.solve()
        acados_u = np.array([acados_solver.get(i, "u") for i in range(qp.N)])
        acados_x = np.array([acados_solver.get(i, "x") for i in range(qp.N+1)])
        acados_lam = np.concatenate([acados_solver.get(i, "lam") for i in range(qp.N+1)])
        acados_pi = np.concatenate([acados_solver.get(i, "pi") for i in range(qp.N)])

        casadi_solver = AcadosCasadiOcpQpSolver(qp)
        status = casadi_solver.solve()
        casadi_u = np.array([casadi_solver.get(i, "u") for i in range(qp.N)])
        casadi_x = np.array([casadi_solver.get(i, "x") for i in range(qp.N+1)])
        casadi_lam = np.concatenate([casadi_solver.get(i, "lam") for i in range(qp.N+1)])
        casadi_pi = np.concatenate([casadi_solver.get(i, "pi") for i in range(qp.N)])

        # evaluate difference
        assert np.allclose(casadi_u, acados_u, atol=1e-5), f"u mismatch for {qp_json_file}"
        assert np.allclose(casadi_x, acados_x, atol=1e-5), f"x mismatch for {qp_json_file}"
        assert np.allclose(casadi_lam, acados_lam, atol=1e-5), f"lam mismatch for {qp_json_file}"
        assert np.allclose(casadi_pi, acados_pi, atol=1e-5), f"pi mismatch for {qp_json_file}"
        print(f"QP solution matches for {qp_json_file}.")

if __name__ == "__main__":
    main()
import os
from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
json_dir = os.path.abspath(os.path.join(script_dir, 'qp_tests'))
qp_json_files = ['prob_0.json',
                 'pendulum_qp.json',
                 'pendulum_slack.json',]
def main():
    for qp_json_file in qp_json_files:
        qp_json_file = os.path.join(json_dir, qp_json_file)
        qp = AcadosOcpQp.from_json(qp_json_file)
        opts = AcadosOcpQpOptions()
        opts.iter_max = 500
        opts.qp_solver = 'PARTIAL_CONDENSING_HPIPM'

        acados_solver = AcadosOcpQpSolver(qp, opts=opts)
        status = acados_solver.solve()
        acados_u = np.array([acados_solver.get(i, "u") for i in range(qp.N)])
        acados_x = np.array([acados_solver.get(i, "x") for i in range(qp.N+1)])
        acados_sl = np.concatenate([acados_solver.get(i, "sl") for i in range(qp.N+1)])
        acados_su = np.concatenate([acados_solver.get(i, "su") for i in range(qp.N+1)])
        acados_lam = np.concatenate([acados_solver.get(i, "lam") for i in range(qp.N+1)])
        acados_pi = np.concatenate([acados_solver.get(i, "pi") for i in range(qp.N)])
        iterate_acados = acados_solver.get_iterate()
        # acados_cost = acados_solver.get_cost()

        casadi_solver = AcadosCasadiOcpQpSolver(qp)
        casadi_solver.set_iterate(iterate_acados) # set initial guess from acados solution
        status = casadi_solver.solve()
        casadi_u = np.array([casadi_solver.get(i, "u") for i in range(qp.N)])
        casadi_x = np.array([casadi_solver.get(i, "x") for i in range(qp.N+1)])
        casadi_sl = np.concatenate([casadi_solver.get(i, "sl") for i in range(qp.N+1)])
        casadi_su = np.concatenate([casadi_solver.get(i, "su") for i in range(qp.N+1)])
        casadi_lam = np.concatenate([casadi_solver.get(i, "lam") for i in range(qp.N+1)])
        casadi_pi = np.concatenate([casadi_solver.get(i, "pi") for i in range(qp.N)])
        # casadi_cost = casadi_solver.get_cost()

        # evaluate difference
        assert np.allclose(casadi_u, acados_u, atol=5e-5, rtol=5e-5), f"u mismatch for {qp_json_file} with error {np.max(np.abs(casadi_u - acados_u))}"
        assert np.allclose(casadi_x, acados_x, atol=5e-5, rtol=5e-5), f"x mismatch for {qp_json_file} with error {np.max(np.abs(casadi_x - acados_x))}"
        assert np.allclose(casadi_sl, acados_sl, atol=5e-5, rtol=5e-5), f"sl mismatch for {qp_json_file} with error {np.max(np.abs(casadi_sl - acados_sl))}"
        assert np.allclose(casadi_su, acados_su, atol=5e-5, rtol=5e-5), f"su mismatch for {qp_json_file} with error {np.max(np.abs(casadi_su - acados_su))}"
        assert np.allclose(casadi_lam, acados_lam, atol=5e-5, rtol=5e-5), f"lam mismatch for {qp_json_file} with error {np.max(np.abs(casadi_lam - acados_lam))}"
        assert np.allclose(casadi_pi, acados_pi, atol=5e-5, rtol=5e-5), f"pi mismatch for {qp_json_file} with error {np.max(np.abs(casadi_pi - acados_pi))}"
        # assert np.isclose(casadi_cost, acados_cost, atol=1e-4), f"cost mismatch for {qp_json_file}"
        print(f"QP solution matches for {qp_json_file}.")

if __name__ == "__main__":
    main()
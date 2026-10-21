import os
from acados_template import AcadosOcpQpSolver, AcadosCasadiOcpQpSolver, AcadosOcpQp, AcadosOcpQpOptions, AcadosOcpQpHpipmOptions
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
json_dir = os.path.abspath(os.path.join(script_dir, 'qp_tests'))
qp_json_files = ['pendulum_qp.json',
                 'pendulum_slack.json',
                 'pend_idxs_rev_min_qp0.json',]
TOL = 5e-5


def get_solution(solver, N: int) -> dict:
    """Get the solution of an OCP QP solver as one stacked array per field."""
    return {
        'u': np.array([solver.get(i, 'u') for i in range(N)]),
        'x': np.array([solver.get(i, 'x') for i in range(N+1)]),
        'sl': np.concatenate([solver.get(i, 'sl') for i in range(N+1)]),
        'su': np.concatenate([solver.get(i, 'su') for i in range(N+1)]),
        'lam': np.concatenate([solver.get(i, 'lam') for i in range(N+1)]),
        'pi': np.concatenate([solver.get(i, 'pi') for i in range(N)]),
    }


def assert_solutions_match(sol: dict, sol_ref: dict, label: str):
    for field, ref in sol_ref.items():
        if sol[field].size == 0 and ref.size == 0:  # no slacks
            continue
        np.testing.assert_allclose(sol[field], ref, atol=TOL, rtol=TOL, err_msg=f"{field} mismatch for {label}")
        print(f"diff in {field} matches for {label} with error {np.max(np.abs(sol[field] - ref)):.2e}")


def solve_with_acados(qp: AcadosOcpQp, opts: AcadosOcpQpOptions) -> AcadosOcpQpSolver:
    solver = AcadosOcpQpSolver(qp, opts=opts)
    status = solver.solve()
    assert status == 0, f"acados QP solver returned status {status}."
    return solver


def check_hpipm_options(qp: AcadosOcpQp, casadi_sol: dict, n_iter_default: int, label: str):
    """Check that the options set in AcadosOcpQpHpipmOptions are used by HPIPM."""
    # same solution, but different number of iterations
    opts = AcadosOcpQpHpipmOptions()
    opts.tol_stat = 1e-6
    opts.tol_eq = 1e-8
    opts.tol_ineq = 1e-8
    opts.tol_comp = 1e-12
    solver = solve_with_acados(qp, opts)
    assert_solutions_match(get_solution(solver, qp.N), casadi_sol, f"{label} (acados, HPIPM options)")
    n_iter = solver.get_stats('iter')
    assert n_iter != n_iter_default, f"HPIPM options did not change the number of iterations ({n_iter}) for {label}."

    # exit condition on minimum step length
    opts = AcadosOcpQpHpipmOptions()
    opts.alpha_min = 0.999
    status = AcadosOcpQpSolver(qp, opts=opts).solve()
    assert status == 3, f"Expected ACADOS_MINSTEP for alpha_min = 0.999 for {label}, got status {status}."


def main():
    for qp_json_file in qp_json_files:
        qp = AcadosOcpQp.from_json(os.path.join(json_dir, qp_json_file))

        # acados, shared options only
        opts = AcadosOcpQpOptions()
        opts.iter_max = 500
        opts.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
        acados_solver = solve_with_acados(qp, opts)

        # CasADi reference, initialized with the acados solution
        casadi_solver = AcadosCasadiOcpQpSolver(qp)
        casadi_solver.set_iterate(acados_solver.get_iterate())
        status = casadi_solver.solve()
        assert status == 0, f"CasADi QP solver returned status {status}."
        casadi_sol = get_solution(casadi_solver, qp.N)

        assert_solutions_match(get_solution(acados_solver, qp.N), casadi_sol, f"{qp_json_file} (acados)")
        check_hpipm_options(qp, casadi_sol, acados_solver.get_stats('iter'), qp_json_file)

if __name__ == "__main__":
    main()

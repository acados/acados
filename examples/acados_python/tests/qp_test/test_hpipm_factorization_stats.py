import pytest

import casadi as ca
import numpy as np

from acados_template import (
    AcadosModel,
    AcadosOcp,
    AcadosOcpQp,
    AcadosOcpQpHpipmOptions,
    AcadosOcpQpSolver,
    AcadosOcpSolver,
)


def _make_box_qp(hessian, linear_term):
    qp = AcadosOcpQp(N=0)
    qp.set("Q", 0, hessian)
    qp.set("q", 0, linear_term)
    qp.set("idxb", 0, np.array([0, 1]))
    qp.set("lbx", 0, np.array([-1.0, -1.0]))
    qp.set("ubx", 0, np.array([1.0, 1.0]))
    qp.set("lbx_mask", 0, np.ones(2))
    qp.set("ubx_mask", 0, np.ones(2))
    qp.set("idxs_rev", 0, np.array([-1, -1]))
    return qp


def _make_ill_conditioned_diagonal_box_qp(objective_scale=1.0):
    hessian = objective_scale * np.diag([1.0, 1e-12])
    linear_term = objective_scale * np.array([-2.0, -1e-8])
    return _make_box_qp(hessian, linear_term)


def _nondiagonal_cost_data():
    angle = 1.0
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    hessian = rotation @ np.diag([1.0, 1e-8]) @ rotation.T
    objective_scale = 1e12
    target = np.array([2.0, 0.25])
    return objective_scale * hessian, target


def _make_ill_conditioned_nondiagonal_box_qp():
    hessian, target = _nondiagonal_cost_data()
    return _make_box_qp(
        hessian,
        -hessian @ target,
    )


def _print_hpipm_diagnostics(solver):
    solver.print_statistics()
    stats = solver.get_stats("statistics")
    iterations = solver.get_stats("iter")
    iteration_stats = stats[1:iterations + 1]
    lq_fact = iteration_stats[:, 13]
    itref_corr = iteration_stats[:, 15]
    linear_residuals = iteration_stats[:, 16:20]

    print("\nHPIPM factorization/refinement diagnostics")
    print(f"LQ iterations: {np.count_nonzero(lq_fact == 1.0)} / {iterations}")
    print(f"Corrector refinement steps: {int(np.sum(itref_corr))}")
    print("iter  lq_fact  itref_corr  lin_res_stat  lin_res_eq  lin_res_ineq  lin_res_comp")
    for iteration, (lq_used, refinement_count, residuals) in enumerate(
        zip(lq_fact, itref_corr, linear_residuals), start=1
    ):
        print(
            f"{iteration:4d}  {int(lq_used):7d}  {int(refinement_count):10d}  "
            f"{residuals[0]:12.3e}  {residuals[1]:10.3e}  "
            f"{residuals[2]:12.3e}  {residuals[3]:12.3e}"
        )

    return iteration_stats


def _make_hpipm_solver(qp, hpipm_mode):
    opts = AcadosOcpQpHpipmOptions()
    opts.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    opts.hpipm_mode = hpipm_mode
    opts.print_level = 1
    return AcadosOcpQpSolver(qp, opts=opts)


def test_robust_mode_reports_lq_factorization():
    print("\n=== Demonstration: ROBUST mode reports LQ factorization ===")
    solver = _make_hpipm_solver(_make_ill_conditioned_diagonal_box_qp(), "ROBUST")
    status = solver.solve()
    assert status == 0

    iteration_stats = _print_hpipm_diagnostics(solver)
    lq_fact = iteration_stats[:, 13]
    itref_corr = iteration_stats[:, 15]
    assert np.any(lq_fact == 1.0), (
        "ROBUST mode should report LQ factorization in statistics column 13; "
        f"lq_fact={lq_fact}, itref_corr={itref_corr}"
    )
    assert np.all(itref_corr >= 0.0)


def test_balance_mode_reports_corrector_refinement_for_scaled_diagonal_qp():
    print("\n=== Demonstration: BALANCE mode reports corrector refinement for a scaled diagonal QP ===")
    solver = _make_hpipm_solver(_make_ill_conditioned_diagonal_box_qp(objective_scale=1e10), "BALANCE")
    status = solver.solve()
    assert status == 0
    iteration_stats = _print_hpipm_diagnostics(solver)
    lq_fact = iteration_stats[:, 13]
    itref_corr = iteration_stats[:, 15]
    assert np.any(itref_corr > 0.0), (
        "BALANCE mode should report at least one corrector refinement; "
        f"itref_corr={itref_corr}, lq_fact={lq_fact}"
    )
    assert np.all(itref_corr >= 0.0)


def test_balance_mode_reports_lq_fallback_for_ill_conditioned_nondiagonal_qp():
    print("\n=== Demonstration: BALANCE mode reports LQ fallback for an ill-conditioned non-diagonal QP ===")
    solver = _make_hpipm_solver(_make_ill_conditioned_nondiagonal_box_qp(), "BALANCE")
    status = solver.solve()
    assert status == 0

    iteration_stats = _print_hpipm_diagnostics(solver)
    lq_fact = iteration_stats[:, 13]
    assert np.any(lq_fact == 1.0), (
        "BALANCE mode should report at least one residual-triggered LQ fallback; "
        f"lq_fact={lq_fact}"
    )


@pytest.mark.parametrize("lq_fact", [-1, 0, 1, 2])
def test_hpipm_sqp_lq_fallback(lq_fact, tmp_path):
    print(f"\n=== Demonstration: terminal-only SQP with lq_fact={lq_fact} ===")
    hessian, target = _nondiagonal_cost_data()

    model = AcadosModel()
    model.name = "balance_lq_terminal_sqp"
    model.x = ca.SX.sym("x", 2)

    ocp = AcadosOcp()
    ocp.name = model.name
    ocp.model = model
    ocp.solver_options.N_horizon = 0
    ocp.solver_options.tf = 1.0
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.nlp_solver_max_iter = 10
    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.hpipm_mode = "BALANCE"
    ocp.solver_options.qp_solver_itref_corr_max = 0
    ocp.solver_options.qp_solver_lq_fact = lq_fact
    ocp.solver_options.qp_tol = 1e-10
    ocp.solver_options.qp_solver_iter_max = 100
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.print_level = 3
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.model.cost_expr_ext_cost_e = 0.5 * (model.x - target).T @ hessian @ (model.x - target)
    ocp.constraints.idxbx_e = np.array([0, 1])
    ocp.constraints.lbx_e = np.array([-1.0, -1.0])
    ocp.constraints.ubx_e = np.array([1.0, 1.0])

    ocp.code_gen_options.code_export_directory = str(tmp_path / "generated_code")
    solver = AcadosOcpSolver(
        ocp,
        verbose=False,
    )
    status = solver.solve()
    solver.print_statistics()
    assert status == 0, f"SQP solver returned status {status}."

    lq_iterations = solver.get_stats("lq")
    corrector_refinements = solver.get_stats("itref_c")
    qp_iterations = solver.get_stats("qp_iter")
    lq_sum = solver.get_stats("hpipm_lq_sum")
    itref_c_sum = solver.get_stats("hpipm_itref_c_sum")
    print(f"SQP iterations: {solver.get_stats('nlp_iter')}")
    print(f"LQ factorizations per SQP iteration: {lq_iterations}")
    print(f"Corrector refinement steps per SQP iteration: {corrector_refinements}")
    print(f"QP iterations per SQP iteration: {qp_iterations}")
    assert np.isclose(lq_sum, np.sum(lq_iterations[1:])), (
        f"hpipm_lq_sum={lq_sum} does not match the sum of per-iteration values "
        f"{np.sum(lq_iterations[1:])}"
    )
    assert np.isclose(itref_c_sum, np.sum(corrector_refinements[1:])), (
        f"hpipm_itref_c_sum={itref_c_sum} does not match the sum of per-iteration values "
        f"{np.sum(corrector_refinements[1:])}"
    )
    assert np.all(corrector_refinements[1:] == 0), (
        "qp_solver_itref_corr_max=0 should disable corrector refinement; "
        f"itref_c={corrector_refinements}"
    )
    if lq_fact in (-1, 1):
        assert np.any(lq_iterations[1:] > 0), (
            f"lq_fact={lq_fact} with BALANCE mode should retain residual-triggered LQ fallback; "
            f"lq={lq_iterations}"
        )
    elif lq_fact == 0:
        assert np.all(lq_iterations[1:] == 0), (
            "lq_fact=0 should disable LQ factorization; "
            f"lq={lq_iterations}"
        )
    else:
        assert np.all(lq_iterations[1:] == qp_iterations[1:]), (
            "lq_fact=2 should use LQ for every HPIPM iteration; "
            f"lq={lq_iterations}, qp_iter={qp_iterations}"
        )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))

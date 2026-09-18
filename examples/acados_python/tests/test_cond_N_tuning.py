# Test for cond_N_tuning: measured choice of qp_solver_cond_N.
# Run from examples/acados_python/tests, like the other tests here:
#     python test_cond_N_tuning.py
#
# What is checked is BEHAVIOUR, not speed (the CI machine's timings are its own):
#   * the offline tuner returns a horizon from the candidate set, sets it on the
#     solver, and leaves the solver's iterate and solution unchanged;
#   * the solution at the chosen horizon equals the solution at the default
#     (condensing changes the cost of the QP, never its solution);
#   * with a tolerance of 100 % nothing is changed;
#   * the online tuner commits within one probe per candidate and the closed
#     loop it drives produces the same controls as the default;
#   * a full-condensing solver is refused with a clear error.
import sys
sys.path.insert(0, '../pendulum_on_cart/common')
import numpy as np
import scipy.linalg
from acados_template import AcadosOcp, AcadosOcpSolver
from acados_template.cond_N_tuning import tune_qp_solver_cond_N, CondNTuner, candidates
from pendulum_model import export_pendulum_ode_model


def make_solver(qp_solver='PARTIAL_CONDENSING_HPIPM', N=100, name='cond_N_tuning_test'):
    ocp = AcadosOcp()
    model = export_pendulum_ode_model()
    model.name = name
    ocp.model = model
    nx, nu = model.x.rows(), model.u.rows()
    ny = nx+nu
    ocp.solver_options.N_horizon = N
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])
    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'
    ocp.cost.Vx = np.zeros((ny, nx)); ocp.cost.Vx[:nx, :nx] = np.eye(nx)
    ocp.cost.Vu = np.zeros((ny, nu)); ocp.cost.Vu[nx:, :] = np.eye(nu)
    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.constraints.lbu = np.array([-80.])
    ocp.constraints.ubu = np.array([80.])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.nlp_solver_max_iter = 200
    # The same swing-up as test_qp_solver.py (tf = 1), only discretised finer
    # so the horizon is long enough for condensing to matter. Stretching tf
    # instead makes the full-step SQP cycle and hit max_iter.
    ocp.solver_options.tf = 1.0
    return AcadosOcpSolver(ocp, json_file=f'{name}.json'), N


def cold_start(solver, x0):
    """Iterate reset with every stage at x0 (from zero the full-step SQP
    diverges on this problem). The reset takes x0 from lbx_0, so set it first."""
    solver.set(0, 'lbx', x0)
    solver.set(0, 'ubx', x0)
    solver.reset(reset_x_to_x0_bar=True)


def test_offline():
    solver, N = make_solver()
    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    solver.solve_for_x0(x0_bar=x0)
    before = solver.store_iterate_to_flat_obj()
    default = int(solver.ocp.solver_options.qp_solver_cond_N)
    assert default == N, f'expected acados default cond_N = N, got {default}'

    choice = tune_qp_solver_cond_N(solver, verbose=True)
    assert choice in candidates(N), f'choice {choice} not a candidate'
    assert int(solver.ocp.solver_options.qp_solver_cond_N) == choice, 'cond_N not set on the solver'
    after = solver.store_iterate_to_flat_obj()
    assert before.allclose(after, rtol=0., atol=1e-12), 'tuning changed the solver iterate'

    # same solution from the same start at the chosen horizon and at the default
    # (both SQP runs stop at nlp tol 1e-6, so equal up to that)
    cold_start(solver, x0); solver.solve_for_x0(x0_bar=x0); tuned = solver.store_iterate_to_flat_obj()
    solver.update_qp_solver_cond_N(N)
    cold_start(solver, x0); solver.solve_for_x0(x0_bar=x0); reference = solver.store_iterate_to_flat_obj()
    assert reference.allclose(tuned, atol=1e-5), f'solution at cond_N={choice} differs from default'

    # nothing is better by more than 100 %: the current value must be kept
    solver.update_qp_solver_cond_N(N)
    solver.solve_for_x0(x0_bar=x0)
    kept = tune_qp_solver_cond_N(solver, tolerance=1.0)
    assert kept == N, f'with tolerance 1.0 the tuner must keep cond_N = N, got {kept}'
    print(f'test_offline ok: chose cond_N={choice} (default {N})')


def test_online():
    solver, N = make_solver(name='cond_N_tuning_test_online')
    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    steps = 3*len(candidates(N))
    # default closed loop (open-loop plant: apply the first control through the OCP's own model)
    def loop(with_tuner):
        tuner = CondNTuner(solver) if with_tuner else None
        x, us = x0.copy(), []
        for _ in range(steps):
            if tuner: tuner.before_solve()
            u = solver.solve_for_x0(x0_bar=x)
            if tuner: tuner.after_solve()
            us.append(u.copy())
            x = solver.get(1, 'x')          # next state from the OCP's own prediction
        return np.array(us), tuner

    solver.update_qp_solver_cond_N(N); cold_start(solver, x0)
    reference, _ = loop(False)
    solver.update_qp_solver_cond_N(N); cold_start(solver, x0)
    tuned, tuner = loop(True)
    assert tuner.done, 'online tuner did not commit'
    assert tuner.choice in candidates(N), tuner.choice
    assert tuner.probes <= len(candidates(N)), f'{tuner.probes} probes for {len(candidates(N))} candidates'
    err = np.abs(tuned-reference).max()/max(np.abs(reference).max(), 1.)
    assert err < 1e-4, f'controls with the online tuner differ from the default by {err:.1e}'
    print(f'test_online ok: committed cond_N={tuner.choice} after {tuner.probes} probes')


def test_refuses_full_condensing():
    solver, N = make_solver('FULL_CONDENSING_HPIPM', N=20, name='cond_N_tuning_test_full')
    solver.solve()
    try:
        tune_qp_solver_cond_N(solver)
    except ValueError as e:
        print(f'test_refuses_full_condensing ok: {e}')
        return
    raise AssertionError('a full-condensing solver must be refused')


if __name__ == '__main__':
    test_offline()
    test_online()
    test_refuses_full_condensing()
    print('all cond_N_tuning tests passed')

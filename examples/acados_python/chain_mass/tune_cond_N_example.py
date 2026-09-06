"""Example: choose `qp_solver_cond_N` by measurement, on the chain-of-masses OCP.

The partial-condensing horizon `qp_solver_cond_N` is set by hand and defaults
to N (no condensing). Its best value depends on N, nx, nu and the machine; on
this example the default costs 1.7-4.3x in QP time. This example builds the
chain OCP as in main.py, solves once, calls `tune_qp_solver_cond_N`, and shows
the closed-loop step time before and after.

    python3 tune_cond_N_example.py            # n_mass=8, N=100
    python3 tune_cond_N_example.py 12 40
"""
import sys
import numpy as np
import scipy.linalg
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSim, AcadosSimSolver
from export_chain_mass_model import export_chain_mass_model
from utils import compute_steady_state
from acados_template.cond_N_tuning import tune_qp_solver_cond_N


def build(n_mass, N, Ts=0.2, m=0.033, D=1.0, L=0.033):
    M = n_mass-2
    model = export_chain_mass_model(n_mass, m, D, L)
    model.name = f'chain_tune_{n_mass}_{N}'
    ocp = AcadosOcp()
    ocp.model = model
    nx, nu = model.x.rows(), model.u.rows()
    ny = nx+nu
    xEndRef = np.zeros((3, 1)); xEndRef[0] = L*(M+1)*6
    xrest = compute_steady_state(n_mass, m, D, L, np.zeros((3, 1)), xEndRef)
    ocp.solver_options.N_horizon = N
    ocp.cost.cost_type, ocp.cost.cost_type_e = 'LINEAR_LS', 'LINEAR_LS'
    q_diag = np.ones((nx, 1)); q_diag[3*M:3*M+3] = M+1
    Q, R = 2*np.diagflat(q_diag), 2*np.diagflat(1e-2*np.ones((nu, 1)))
    ocp.cost.W, ocp.cost.W_e = scipy.linalg.block_diag(Q, R), Q
    ocp.cost.Vx = np.zeros((ny, nx)); ocp.cost.Vx[:nx, :nx] = np.eye(nx)
    ocp.cost.Vu = np.zeros((ny, nu)); ocp.cost.Vu[nx:, :] = np.eye(nu)
    ocp.cost.Vx_e = np.eye(nx)
    ocp.cost.yref = np.vstack((xrest, np.zeros((nu, 1)))).flatten()
    ocp.cost.yref_e = xrest.flatten()
    ocp.constraints.lbu, ocp.constraints.ubu = -np.ones(nu), np.ones(nu)
    ocp.constraints.idxbu = np.arange(nu)
    ocp.constraints.x0 = xrest.reshape((nx,))
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.sim_method_num_stages = 2
    ocp.solver_options.sim_method_num_steps = 2
    ocp.solver_options.qp_tol = 1e-6
    ocp.solver_options.tol = 1e-6
    ocp.solver_options.tf = N*Ts
    solver = AcadosOcpSolver(ocp, json_file=f'{model.name}.json')
    sim = AcadosSim()
    sim.model = export_chain_mass_model(n_mass, m, D, L)
    sim.model.name = f'chain_tune_sim_{n_mass}'
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages, sim.solver_options.num_steps = 2, 2
    sim.solver_options.T = Ts
    integrator = AcadosSimSolver(sim, json_file=f'{sim.model.name}.json')
    x0 = xrest.reshape((nx,))
    for _ in range(5):
        x0 = integrator.simulate(x=x0, u=np.array([-1., 1., 1.]))
    return solver, integrator, x0


def closed_loop(solver, integrator, x0, steps):
    x, tot, tqp = x0.copy(), 0., 0.
    for _ in range(steps):
        u = solver.solve_for_x0(x0_bar=x)
        tot += solver.get_stats('time_tot'); tqp += solver.get_stats('time_qp')
        x = integrator.simulate(x=x, u=u)
    return tot, tqp


if __name__ == '__main__':
    n_mass = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    solver, integrator, x0 = build(n_mass, N)
    steps = 25
    default_tot, default_qp = closed_loop(solver, integrator, x0, steps)     # cond_N = N
    solver.solve_for_x0(x0_bar=x0)
    choice = tune_qp_solver_cond_N(solver, verbose=True)
    tuned_tot, tuned_qp = closed_loop(solver, integrator, x0, steps)
    print(f'\nchain n_mass={n_mass} (nx={solver.get(0, "x").size}), N={N}, {steps} closed-loop steps')
    print(f'  qp_solver_cond_N = {N} (default): {default_tot*1e3:8.1f} ms total, {default_qp*1e3:8.1f} ms QP')
    print(f'  qp_solver_cond_N = {choice} (tuned):  {tuned_tot*1e3:8.1f} ms total, {tuned_qp*1e3:8.1f} ms QP')
    print(f'  -> {default_tot/tuned_tot:.2f}x faster overall, {default_qp/tuned_qp:.2f}x on the QP')

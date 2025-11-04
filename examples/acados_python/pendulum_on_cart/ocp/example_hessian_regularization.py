import sys
sys.path.insert(0, '../common')
import casadi as ca
import numpy as np
from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model

def formulate_ocp(Tf: float = 1.0, N: int = 20, regularize_method: str = 'NO_REGULARIZE') -> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = np.array([0, np.pi, 0, 0])
    ocp.constraints.remove_x0_elimination()

    # set options
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N
    ocp.solver_options.store_iterates = True
    ocp.solver_options.nlp_solver_max_iter = 300
    ocp.solver_options.nlp_solver_tol_stat = 1e-7
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES, PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_HPIPM
    ocp.solver_options.regularize_method = regularize_method
    ocp.solver_options.qp_solver_ric_alg  = 0
    ocp.solver_options.reg_epsilon = 1e-4
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    return ocp

def main(regularize_method):
    print(f"Testing regularization method: {regularize_method}")
    N = 100
    dt = 0.01
    Tf = N * dt
    ocp = formulate_ocp(Tf=Tf, N=N, regularize_method=regularize_method)
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    initial_guess = ocp.create_default_initial_iterate()
    angle_init = np.linspace(np.pi, 0, N)
    x_traj = np.zeros((N+1, 4))
    for i in range(N):
        x_traj[i][1] = angle_init[i]
    initial_guess.x_traj = x_traj
    ocp_solver.load_iterate_from_obj(initial_guess)

    ocp_solver.solve()
    ocp_solver.print_statistics()

    eigs_full = ocp_solver.qp_diagnostics('FULL_HESSIAN')
    eigs_proj = ocp_solver.qp_diagnostics('PROJECTED_HESSIAN')
    print(f"Full Hessian min eigenvalue: {eigs_full['min_eigv_global']}")
    print(f"Projected Hessian min eigenvalue: {eigs_proj['min_eigv_global']}")

    if regularize_method != 'NO_REGULARIZE':
        assert eigs_full['min_eigv_global'] >= 0, "Full Hessian is indefinite!"
        assert eigs_proj['min_eigv_global'] >= 0, "Projected Hessian is indefinite!"
    ocp_solver = None

if __name__ == "__main__":
    main(regularize_method='NO_REGULARIZE')
    main(regularize_method='PROJECT')
    main(regularize_method='MIRROR')
    main(regularize_method='CONVEXIFY')
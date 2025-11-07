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
    ocp.solver_options.nlp_solver_max_iter = 300
    ocp.solver_options.nlp_solver_tol_stat = 1e-7
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES, PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_HPIPM
    ocp.solver_options.regularize_method = regularize_method
    ocp.solver_options.qp_solver_ric_alg  = 0
    ocp.solver_options.reg_epsilon = 1e-4
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    return ocp

def main(regularize_method):
    print(f"Testing regularization method: {regularize_method}")
    N = 10
    dt = 0.05
    Tf = N * dt
    ocp = formulate_ocp(Tf=Tf, N=N, regularize_method=regularize_method)
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    tol = 1e-6
    for i in range(20):
        status = ocp_solver.solve()
        eigs_full = ocp_solver.qp_diagnostics('FULL_HESSIAN')
        eigs_proj = ocp_solver.qp_diagnostics('PROJECTED_HESSIAN')
        if regularize_method != 'NO_REGULARIZE':
            assert eigs_full['min_eigv_global'] >= 0, "Full Hessian is indefinite!"
            assert eigs_proj['min_eigv_global'] >= 0, "Projected Hessian is indefinite!"
        residuals = ocp_solver.get_residuals(recompute=True)
        if max(residuals) < tol or status != 0:
            break
    ocp_solver = None

if __name__ == "__main__":
    main(regularize_method='NO_REGULARIZE')
    main(regularize_method='PROJECT')
    main(regularize_method='MIRROR')
    main(regularize_method='CONVEXIFY')
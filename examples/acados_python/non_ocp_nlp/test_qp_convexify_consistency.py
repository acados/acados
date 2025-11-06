import sys
sys.path.insert(0, '../getting_started')
import casadi as ca
from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np


def formulate_ocp(N=20, Tf=1.0, Fmax=80, reg='NO_REGULARIZE'):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    # # nonlinear
    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q, R).full()

    # terminal cost
    Q_e = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_e

    # set constraints
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    # set options
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N
    ocp.solver_options.nlp_solver_max_iter = 300
    ocp.solver_options.nlp_solver_tol_stat = 1e-7
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES, PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_HPIPM
    ocp.solver_options.reg_epsilon = 1e-4
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.regularize_method = reg  # 'NO_REGULARIZATION', 'PROJECT', 'CONVEXIFY', 'MIRROR'
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.exact_hess_dyn = True
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'

    return ocp

def main(reg):
    N = 20
    dt = 0.01
    Tf = N*dt
    ocp = formulate_ocp(N=N, Tf=Tf, reg=reg)
    print('\nregularization method:{}'.format(reg))

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    initial_guess = ocp.create_default_initial_iterate()
    ocp_solver.load_iterate_from_obj(initial_guess)
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    residuals = ocp_solver.get_residuals(recompute=True)
    cost = ocp_solver.get_cost()
    print('status = ', status)
    print("cost function value of solution = ", cost)
    print(f"NLP residuals = {residuals}")


if __name__ == "__main__":
    main(reg='CONVEXIFY')
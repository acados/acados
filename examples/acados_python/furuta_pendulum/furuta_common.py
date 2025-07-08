from acados_template import AcadosModel
import casadi as ca
import numpy as np

from acados_template import AcadosOcp, AcadosOcpSolver
import scipy.linalg


def get_furuta_model():

    # Distances
    L1 = 0.1035 # 103.5mm
    l2 = 0.0955 # 92.1mm

    # mass
    m2 = 0.192 # 199g
    J2 = 7.653e-04  # kg/mm^2
    g = 9.81  # N/kg

    # inertia arm 1
    J1_ges = 5.3875e-04 + 0.75e-04  # J1 + m1*l1^2

    # inertia arm 2
    J2_ges = J2 + m2*l2**2

    # total inertia at motor
    J0 = J1_ges + m2*L1**2

    # damping hub motor
    b1 = 40*1e-4

    # damping coupling between both arms
    k = 0.098
    b2 = 2*k*J2_ges

    # applied torques
    tau2 = 0

    # named symbolic variables
    theta1 = ca.SX.sym('theta1')  # angle around vertical axis (axis next to motor) [rad]
    theta2 = ca.SX.sym('theta2')  # angle around horizontal axis (axis next to mass) [rad]
    dtheta1 = ca.SX.sym('dtheta1')  # angular velocity of rod 1 [rad/s]
    dtheta2 = ca.SX.sym('dtheta2')  # angular velocity of rod 2 [rad/s]
    dtheta = ca.vertcat(dtheta1, dtheta2)
    tau1 = ca.SX.sym('tau1')  # torque acting on first rod [Nm]

    x = ca.vertcat(theta1, theta2, dtheta1, dtheta2)
    xdot = ca.SX.sym('xdot', x.shape)
    u = tau1
    theta2 = theta2 - np.pi

    # dynamics
    sin_theta_2 = ca.sin(theta2)
    cos_theta_2 = ca.cos(theta2)
    sin_2_theta_2 = ca.sin(2*theta2)

    factor = m2*L1*l2

    Matrix1 = ca.blockcat([
                [J0 + J2_ges*sin_theta_2**2, factor*cos_theta_2],
                [factor*cos_theta_2, J2_ges]])
    Matrix2 = ca.blockcat([
                [b1 + 0.5*dtheta2*J2_ges*sin_2_theta_2, 0.5*dtheta1*J2_ges*sin_2_theta_2 - factor*sin_theta_2*dtheta2],
                [-0.5*dtheta1*J2_ges*sin_2_theta_2, b2]])

    rhs = ca.vertcat(tau1, tau2) - Matrix2 @ ca.vertcat(dtheta1, dtheta2) - ca.vertcat(0, g*m2*l2*sin_theta_2)

    f_expl_expr = ca.vertcat(dtheta, ca.solve(Matrix1, rhs))
    f_impl_expr = ca.vertcat(
        dtheta - xdot[:2],
        Matrix1 @ xdot[2:] - rhs
    )

    model = AcadosModel()
    model.name = 'furuta_model'
    model.x = x
    model.xdot = xdot
    model.u = u
    model.f_impl_expr = f_impl_expr
    model.f_expl_expr = f_expl_expr
    return model


def setup_ocp_solver(x0, umax, dt_0, N_horizon, Tf, RTI=False, timeout_max_time=0.0, heuristic="ZERO", with_anderson_acceleration=False, nlp_solver_max_iter = 20, tol = 1e-6, with_abs_cost=False):
    ocp = AcadosOcp()

    model = get_furuta_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    ocp.solver_options.N_horizon = N_horizon

    # set cost module
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    Q_mat = np.diag([50., 500., 1., 1.])
    R_mat = np.diag([1e3])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    ocp.constraints.lbu = np.array([-umax])
    ocp.constraints.ubu = np.array([+umax])
    ocp.constraints.idxbu = np.array([0])

    if with_abs_cost:
        val = 1.4
        # add cost term abs(x[0]-val) via slacks
        # ocp.constraints.idxbx_e = np.array([0])
        # ocp.constraints.lbx_e = np.array([val])
        # ocp.constraints.ubx_e = np.array([val])
        # ocp.constraints.idxsbx_e = np.array([0])

        # ocp.cost.zl_e = 1e2 * np.array([1.0])
        # ocp.cost.zu_e = 1e2 * np.array([1.0])
        # ocp.cost.Zl_e = np.array([10.0])
        # ocp.cost.Zu_e = np.array([10.0])

        ocp.constraints.idxbx = np.array([0])
        ocp.constraints.lbx = np.array([val])
        ocp.constraints.ubx = np.array([val])
        ocp.constraints.idxsbx = np.array([0])

        ocp.cost.zl = 1e3 * np.array([1.0])
        ocp.cost.zu = 1e3 * np.array([1.0])
        ocp.cost.Zl = 0.0 * np.array([1.0])
        ocp.cost.Zu = 0.0 * np.array([1.0])


    ocp.constraints.x0 = x0

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'

    # NOTE we use a nonuniform grid!
    ocp.solver_options.time_steps = np.array([dt_0] + [(Tf-dt_0)/(N_horizon-1)]*(N_horizon-1))
    ocp.solver_options.sim_method_num_steps = np.array([1] + [2]*(N_horizon-1))
    ocp.solver_options.levenberg_marquardt = 1e-6
    ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter
    ocp.solver_options.with_anderson_acceleration = with_anderson_acceleration

    ocp.solver_options.nlp_solver_type = 'SQP_RTI' if RTI else 'SQP'
    ocp.solver_options.qp_solver_cond_N = N_horizon
    ocp.solver_options.tol = tol

    ocp.solver_options.tf = Tf

    # timeout
    ocp.solver_options.timeout_max_time = timeout_max_time
    ocp.solver_options.timeout_heuristic = heuristic

    solver_json = 'acados_ocp_' + model.name + '.json'
    ocp_solver = AcadosOcpSolver(ocp, json_file = solver_json, verbose=False)

    return ocp_solver


if __name__ == "__main__":
    get_furuta_model()

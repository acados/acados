from acados_template import AcadosModel
import casadi as ca
import numpy as np

from scipy.linalg import block_diag

def get_furuta_model():

    # system dimensions
    nx = 4
    nu = 1

    # system parameters

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

    # (unnamed) symbolic variables
    x = ca.vertcat(theta1, theta2, dtheta1, dtheta2)
    xdot = ca.SX.sym('xdot', nx, 1)
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

    # # constraints
    # expr_h = u

    # # cost
    # W_x = np.diag([50, 500, 1, 1])
    # W_u = 1000
    # expr_ext_cost_e = x.T @ W_x @ x
    # expr_ext_cost = expr_ext_cost_e + u.T @ W_u @ u

    # nonlinear least squares
    # cost_expr_y = ca.vertcat(x, u)
    # W = np.blo(W_x, W_u)
    # model.cost_expr_y_e = x
    # model.W_e = W_x * 10

    # populate structure
    model = AcadosModel()
    model.name = 'furuta_model'
    model.x = x
    model.xdot = xdot
    model.u = u
    model.f_impl_expr = f_impl_expr
    model.f_expl_expr = f_expl_expr
    return model
    # model.expr_h = expr_h
    # model.expr_ext_cost = expr_ext_cost
    # model.expr_ext_cost_e = expr_ext_cost_e

    # model.cost_expr_y = cost_expr_y
    # model.W = W


    # LQR
    # fd = sqrt(eps)
    # x0 = Matrix([0, 0, 0, 0])
    # u0 = 0
    # A1 = (f(x0 + Matrix([1, 0, 0, 0])*fd, u0) - f(x0 - Matrix([1, 0, 0, 0])*fd, u0)) / (2*fd)
    # A2 = (f(x0 + Matrix([0, 1, 0, 0])*fd, u0) - f(x0 - Matrix([0, 1, 0, 0])*fd, u0)) / (2*fd)
    # A3 = (f(x0 + Matrix([0, 0, 1, 0])*fd, u0) - f(x0 - Matrix([0, 0, 1, 0])*fd, u0)) / (2*fd)
    # A4 = (f(x0 + Matrix([0, 0, 0, 1])*fd, u0) - f(x0 - Matrix([0, 0, 0, 1])*fd, u0)) / (2*fd)

    # B = (f(x0, u0 + fd) - f(x0, u0 - fd)) / (2*fd)
    # A = Matrix([[A1, A2, A3, A4]])
    # Ts = 0.025
    # K_LQR, P_LQR, eig_of_ctrl_system = lqrd(A, B, W_x, W_u, Ts)

if __name__ == "__main__":
    get_furuta_model()

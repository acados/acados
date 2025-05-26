import casadi as ca
import numpy as np
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, AcadosOcpFlattenedIterate, AcadosCasadiOcpSolver, ACADOS_INFTY


NX = 4
GOAL_POSITION_RADIUS = (5e-02)**2
LENGTH_PENDULUM = 0.8

def create_pendulum_model():
    M = 1.0
    m = 0.1
    g = 9.81
    states = ca.MX.sym("x", NX)
    states_dot = ca.MX.sym('xdot', NX)
    [p, pDot, theta, thetaDot] = ca.vertsplit(states)
    controls = ca.MX.sym("u")
    F = controls
    denominator = M + m - m*ca.cos(theta)*ca.cos(theta)
    f_x_p = pDot
    f_x_pDot = (-m*LENGTH_PENDULUM*ca.sin(theta)*thetaDot*thetaDot + m*g*ca.cos(theta)*ca.sin(theta)+F)/denominator
    f_x_theta = thetaDot
    f_x_thetaDot = (-m*LENGTH_PENDULUM*ca.cos(theta)*ca.sin(theta)*thetaDot*thetaDot + F*ca.cos(theta)+(M+m)*g*ca.sin(theta))/(LENGTH_PENDULUM*denominator)

    f_x = ca.vertcat(f_x_p, f_x_pDot, f_x_theta, f_x_thetaDot)
    f_impl = states_dot - f_x

    model = AcadosModel()
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_x
    model.x = states
    model.xdot = states_dot
    model.u = controls
    model.name = 'pendulum'
    return model

def pendulum_position(p, theta):
    return ca.vertcat(p-LENGTH_PENDULUM*ca.sin(theta), LENGTH_PENDULUM*ca.cos(theta))

def pendulum_final_position_constraint(p, theta):
    pendulum_final_position = pendulum_position(p, theta)
    return ca.sumsqr(pendulum_final_position - ca.vertcat(LENGTH_PENDULUM, LENGTH_PENDULUM)) - GOAL_POSITION_RADIUS

def build_acados_test_problem(mode='GN', with_anderson_acceleration=False, globalization="FIXED_STEP", max_iter=400) -> AcadosOcp:
    print(f"Building acados test problem with mode {mode} and with_anderson_acceleration {with_anderson_acceleration}, {globalization}")

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = create_pendulum_model()

    ocp.model = model

    # Define cost
    R_mat = np.array([[1e-4]])

    ocp.cost.cost_type_0 = 'LINEAR_LS'
    ocp.cost.W_0 = R_mat
    ocp.cost.Vx_0 = np.zeros((1, NX))
    ocp.cost.Vu_0 = np.array([[1]])
    ocp.cost.yref_0 = np.zeros((1, ))

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.W = R_mat
    ocp.cost.Vx = np.zeros((1, NX))
    ocp.cost.Vu = np.array([[1]])
    ocp.cost.yref  = np.zeros((1, ))

    # Define constraints
    ocp.constraints.x0 = np.array([0.0, 0.0, np.pi, 0.0])

    # Convex over Nonlinear Constraints
    if mode in ['GN', 'EXACT']:
        ocp.model.con_h_expr_e = pendulum_final_position_constraint(model.x[0], model.x[2])
        ocp.constraints.lh_e = np.array([-ACADOS_INFTY])
        ocp.constraints.uh_e = np.array([0.0])

    elif mode == 'SCQP':
        r = ca.MX.sym('r', 2, 1)
        ocp.model.con_phi_expr_e = ca.sumsqr(r)
        ocp.model.con_r_in_phi_e = r
        ocp.model.con_r_expr_e = pendulum_position(model.x[0], model.x[2]) - ca.vertcat(LENGTH_PENDULUM, LENGTH_PENDULUM)
        ocp.constraints.lphi_e = np.array([-ACADOS_INFTY])
        ocp.constraints.uphi_e = np.array([GOAL_POSITION_RADIUS])
    else:
        raise ValueError("Wrong mode name!")

    # set solver options
    N_horizon = 20
    dt = 0.05

    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_max_iter = 400
    ocp.solver_options.qp_solver_iter_max = 1000
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.sim_method_num_steps = 5
    ocp.solver_options.sim_method_num_stages = 4
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = dt*N_horizon
    ocp.solver_options.tol = 1e-12
    ocp.solver_options.qp_tol = 1e-1 * ocp.solver_options.tol
    # ocp.solver_options.cost_scaling = np.ones((N_horizon+1, ))
    ocp.solver_options.with_anderson_acceleration = with_anderson_acceleration
    ocp.solver_options.globalization = globalization
    ocp.solver_options.qp_solver_ric_alg = 0
    ocp.solver_options.qp_solver_cond_ric_alg = 0
    ocp.solver_options.nlp_solver_max_iter = max_iter
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_DAQP'
    ocp.solver_options.hessian_approx = 'EXACT' if mode == "EXACT" else 'GAUSS_NEWTON'
    if mode == "EXACT":
        ocp.solver_options.exact_hess_dyn = 1

    if max_iter < 10:
        ocp.solver_options.print_level = 4

    return ocp
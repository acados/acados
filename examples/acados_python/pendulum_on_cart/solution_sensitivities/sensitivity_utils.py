import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, latexify_plot

latexify_plot()

def export_pendulum_ode_model_with_mass_as_param(dt) -> AcadosModel:

    model_name = 'pendulum_parametric'

    # constants
    m = 0.1 # mass of the ball [kg]
    g = 9.81 # gravity constant [m/s^2]
    l = 0.8 # length of the rod [m]

    # set up states & controls
    x1      = ca.SX.sym('x1')
    theta   = ca.SX.sym('theta')
    v1      = ca.SX.sym('v1')
    dtheta  = ca.SX.sym('dtheta')

    x = ca.vertcat(x1, theta, v1, dtheta)

    F = ca.SX.sym('F')
    u = F

    # xdot
    x1_dot      = ca.SX.sym('x1_dot')
    theta_dot   = ca.SX.sym('theta_dot')
    v1_dot      = ca.SX.sym('v1_dot')
    dtheta_dot  = ca.SX.sym('dtheta_dot')

    xdot = ca.vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot)

    # parameters
    m_cart = ca.SX.sym('m_cart')  # mass of the cart [kg]
    p = m_cart

    # dynamics
    cos_theta = ca.cos(theta)
    sin_theta = ca.sin(theta)
    denominator = m_cart + m - m*cos_theta*cos_theta
    f_expl = ca.vertcat(v1,
                       dtheta,
                       (-m*l*sin_theta*dtheta*dtheta + m*g*cos_theta*sin_theta+F)/denominator,
                       (-m*l*cos_theta*sin_theta*dtheta*dtheta + F*cos_theta+(m_cart+m)*g*sin_theta)/(l*denominator)
                       )

    f_impl = xdot - f_expl

    # discrete dynamics via RK4
    ode = ca.Function('ode', [x, u, p], [f_expl])
    k1 = ode(x, u, p)
    k2 = ode(x + dt/2*k1, u, p)
    k3 = ode(x + dt/2*k2, u, p)
    k4 = ode(x + dt*k3, u, p)

    xf = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

    model = AcadosModel()

    model.disc_dyn_expr = xf
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p
    model.name = model_name

    return model


def export_parametric_ocp(
    x0=np.array([0.0, np.pi / 6, 0.0, 0.0]), N_horizon=50, T_horizon=2.0, Fmax=80.0,
    hessian_approx = "GAUSS_NEWTON", qp_solver_ric_alg=1
) -> AcadosOcp:

    ocp = AcadosOcp()
    dt = T_horizon/N_horizon
    ocp.model = export_pendulum_ode_model_with_mass_as_param(dt=dt)

    ocp.dims.N = N_horizon

    Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2 * np.diag([1e-1])

    ocp.cost.cost_type = "EXTERNAL"
    ocp.cost.cost_type_e = "EXTERNAL"

    # NOTE here we make the cost parametric
    ocp.model.cost_expr_ext_cost = ca.exp(ocp.model.p) * ocp.model.x.T @ Q_mat @ ocp.model.x + ocp.model.u.T @ R_mat @ ocp.model.u
    ocp.model.cost_expr_ext_cost_e = ca.exp(ocp.model.p) * ocp.model.x.T @ Q_mat @ ocp.model.x

    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = x0

    # set mass to one
    ocp.parameter_values = np.ones((1,))

    ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.nlp_solver_type = "SQP"
    ocp.solver_options.qp_solver_cond_N = N_horizon

    ocp.solver_options.tf = T_horizon

    ocp.solver_options.qp_solver_ric_alg = qp_solver_ric_alg
    ocp.solver_options.hessian_approx = hessian_approx
    if hessian_approx == 'EXACT':
        ocp.solver_options.nlp_solver_step_length = 0.0
        ocp.solver_options.nlp_solver_max_iter = 1
        ocp.solver_options.qp_solver_iter_max = 200
        ocp.solver_options.tol = 1e-10
        ocp.solver_options.with_solution_sens_wrt_params = True
        ocp.solver_options.with_value_sens_wrt_params = True
    else:
        ocp.solver_options.nlp_solver_max_iter = 400
        ocp.solver_options.tol = 1e-8

    return ocp


def evaluate_hessian_eigenvalues(acados_solver: AcadosOcpSolver, N_horizon: int):
    offset = 0
    min_eigv_total = 1e12
    min_abs_eigv = 1e12

    for i in range(N_horizon+1):
        hess_block_acados = acados_solver.get_hessian_block(i)
        nv = hess_block_acados.shape[0]
        offset += nv

        eigv = np.linalg.eigvals(hess_block_acados)
        min_eigv = np.min(eigv)
        min_eigv_total = min(min_eigv, min_eigv_total)
        min_abs_eigv = min(min_abs_eigv, np.min(np.abs(eigv)))

    # check projected Hessian
    min_abs_eig_proj_hess = 1e12
    min_eig_proj_hess = 1e12
    min_eig_P = 1e12
    min_abs_eig_P = 1e12
    for i in range(1, N_horizon):
        P_mat = acados_solver.get_from_qp_in(i, 'P')
        B_mat = acados_solver.get_from_qp_in(i-1, 'B')
        # Lr: lower triangular decomposition of R within Riccati != R in qp_in!
        Lr = acados_solver.get_from_qp_in(i-1, 'Lr')
        R_ric = Lr @ Lr.T
        proj_hess_block = R_ric + B_mat.T @ P_mat @ B_mat
        eigv = np.linalg.eigvals(proj_hess_block)
        min_eigv = np.min(eigv)
        min_eig_proj_hess = min(min_eigv, min_eig_proj_hess)
        min_abs_eig_proj_hess = min(min_abs_eig_proj_hess, np.min(np.abs(eigv)))
        # P
        eigv = np.linalg.eigvals(P_mat)
        min_eig_P = min(min_eig_P, np.min(eigv))
        min_abs_eig_P = min(min_abs_eig_P, np.min(np.abs(eigv)))

    return min_eigv_total, min_abs_eigv, min_abs_eig_proj_hess, min_eig_proj_hess, min_eig_P, min_abs_eig_P


def plot_cost_gradient_results(p_test, cost_values, acados_cost_grad, np_cost_grad, cost_reconstructed_np_grad):
    _, ax = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(9,9))

    ax[0].plot(p_test, cost_values, label='cost acados', color='k')
    ax[0].plot(p_test, cost_reconstructed_np_grad, "--", label='reconstructed from finite diff')
    ax[0].set_ylabel(r"cost")

    ax[1].plot(p_test, np_cost_grad, "--", label='finite diff')
    ax[1].plot(p_test, acados_cost_grad, "--", label='acados')
    ax[1].set_ylabel(r"$\partial_p V^*$")
    ax[1].set_yscale("log")

    # plot differences
    isub = 2
    ax[isub].plot(p_test, np.abs(acados_cost_grad - np_cost_grad), "--", label='acados vs. finite diff')
    ax[isub].set_ylabel(r"difference $\partial_p V^*$")
    ax[isub].set_yscale("log")

    isub += 1
    ax[isub].plot(p_test, np.abs(acados_cost_grad - np_cost_grad) / np.abs(np_cost_grad), "--", label='acados vs. finite diff')
    ax[isub].set_ylabel(r"rel. diff. $\partial_p V^*$")
    ax[isub].set_yscale("log")

    for i in range(isub+1):
        ax[i].grid()
        ax[i].legend()

    ax[-1].set_xlabel(f"mass")
    ax[-1].set_xlim([p_test[0], p_test[-1]])

    fig_filename = f"cost_gradient.pdf"
    plt.savefig(fig_filename)
    print(f"stored figure as {fig_filename}")
    plt.show()


def plot_results(p_test, pi, pi_reconstructed_acados, pi_reconstructed_np_grad, sens_u, np_grad,
                 min_eig_full=None, min_eig_proj_hess=None, min_eig_P=None,
                 min_abs_eig_full=None, min_abs_eig_proj_hess=None, min_abs_eig_P=None,
                 eigen_analysis=False, qp_solver_ric_alg=1, parameter_name=""):

    nsub = 5 if eigen_analysis else 3

    _, ax = plt.subplots(nrows=nsub, ncols=1, sharex=True, figsize=(9,9))

    isub = 0
    ax[isub].plot(p_test, pi, label='acados', color='k')
    ax[isub].plot(p_test, pi_reconstructed_acados, "--", label='reconstructed from acados')
    ax[isub].plot(p_test, pi_reconstructed_np_grad, "--", label='reconstructed from finite diff')
    ax[isub].set_ylabel(r"$u$")
    ax[isub].set_title(f'qp_solver_ric_alg {qp_solver_ric_alg}')

    isub += 1
    ax[isub].plot(p_test, sens_u, label="acados")
    ax[isub].plot(p_test, np_grad, "--", label="finite diff")
    ax[isub].set_xlim([p_test[0], p_test[-1]])
    ax[isub].set_ylabel(r"$\partial_p u$")

    isub += 1
    ax[isub].plot(p_test, np.abs(sens_u- np_grad), "--", label='acados - finite diff')
    ax[isub].set_ylabel(r"diff $\partial_p u$")
    ax[isub].set_yscale("log")

    if eigen_analysis:
        isub += 1
        ax[isub].plot(p_test, np.sign(min_eig_full), linestyle="-", alpha=.6, label='full Hessian')
        ax[isub].plot(p_test, np.sign(min_eig_proj_hess), "--", label='proj Hessian')
        ax[isub].plot(p_test, np.sign(min_eig_P), ":", label='$P$ Riccati')
        ax[isub].set_ylabel("sign min eig")

        isub += 1
        ax[isub].plot(p_test, min_abs_eig_full, "--", label='full Hessian')
        ax[isub].plot(p_test, min_abs_eig_proj_hess, "--", label='proj Hessian')
        ax[isub].plot(p_test, min_abs_eig_P, "--", label='$P$ Riccati')
        ax[isub].set_ylabel("min abs eig")
        ax[isub].set_yscale("log")

    for i in range(isub+1):
        ax[i].grid()
        ax[i].legend()

    ax[-1].set_xlabel(f"{parameter_name}")

    fig_filename = f"solution_sens_{qp_solver_ric_alg}.pdf"
    plt.savefig(fig_filename)
    print(f"stored figure as {fig_filename}")
    plt.show()

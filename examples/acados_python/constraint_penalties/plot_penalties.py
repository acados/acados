import casadi as ca
import matplotlib.pyplot as plt
import numpy as np

from acados_template import latexify_plot, symmetric_huber_penalty, one_sided_huber_penalty, AcadosOcp

def plot_huber_penalty(symmetric=True):
    u = ca.SX.sym("u")

    delta = 1e-1
    tau = 1e4

    if symmetric:
        penalty, penalty_grad, penalty_hess, penalty_hess_xgn = symmetric_huber_penalty(
            u, delta, tau=tau
        )
    else:
        penalty, penalty_grad, penalty_hess, penalty_hess_xgn = one_sided_huber_penalty(
            u, delta, tau=tau
        )

    huber_penalty_fun = ca.Function("penalty", [u], [penalty])

    u_ = ca.SX.sym("u_")
    huber_penalty_quad_fun = ca.Function(
        "penalty",
        [u_, u],
        [penalty + penalty_grad * (u_ - u) + 0.5 * penalty_hess_xgn * (u_ - u) ** 2],
    )

    us = np.linspace(-1.5, 1.5, 800)
    ys = huber_penalty_fun(us).full()
    latexify_plot()
    plt.figure(figsize=(10, 10))
    plt.plot(us, ys, label="Huber penalty")

    u_lins = [-0.999, -0.8, 0.4]
    for i in range(len(u_lins)):
        ys_quad = huber_penalty_quad_fun(us, u_lins[i]).full()
        plt.plot(
            us, ys_quad, "--", label=f"quadratic approximation {i} at {u_lins[i]:.1f}"
        )

    plt.grid()
    plt.legend()
    plt.xlim(us[0], us[-1])
    plt.ylim(min(ys), max(ys))
    plt.show()


def plot_huber_through_ocp():

    delta = 1e-1
    weight = 1e2
    lower = -10
    upper = 10

    # create dummy ocp
    ocp = AcadosOcp()
    ocp.cost.cost_type = "CONVEX_OVER_NONLINEAR"
    ocp.model.cost_y_expr = ca.SX.zeros(0)
    ocp.model.cost_psi_expr = ca.SX.zeros(1)
    ocp.model.cost_r_in_psi_expr = ca.SX.zeros(0)
    x = ca.SX.sym('x', 1)
    ocp.model.x = x
    ocp.formulate_constraint_as_Huber_penalty(x, weight, upper_bound=upper, lower_bound=lower, huber_delta=delta)

    # extract penalty
    penalty_expr = ca.substitute(ocp.model.cost_psi_expr,
                                 ocp.model.cost_r_in_psi_expr,
                                 ocp.model.cost_y_expr - ocp.cost.yref)
    huber_penalty_fun = ca.Function('penalty', [x], [penalty_expr])

    # extract quadratic approximation
    x_ = ca.SX.sym('x_')
    penalty_hess_xgn = ca.substitute(ocp.model.cost_conl_custom_outer_hess, ocp.model.cost_r_in_psi_expr, ocp.model.cost_y_expr - ocp.cost.yref)
    penalty_hess_ggn, penalty_grad = ca.hessian(penalty_expr, x)
    # x linearization point, x_ for evaluation
    huber_ggn_quad_fun = ca.Function('penalty', [x_, x],
                                     [penalty_expr + penalty_grad*(x_ - x) + 0.5*penalty_hess_ggn*(x_ - x)**2])

    # compute normalized constraint expression
    width = upper - lower
    center = lower + 0.5 * width
    normalized_constr_expr = 2 * (x - center) / width
    J_normalized = ca.jacobian(normalized_constr_expr, x)

    huber_xgn_quad_fun = ca.Function('penalty', [x_, x], [penalty_expr + penalty_grad*(x_ - x) + 0.5*J_normalized.T@penalty_hess_xgn@J_normalized*(x_ - x)**2])

    # create plot data
    lower_plot = 2.5 *lower
    upper_plot = 2.5 *upper
    us = np.linspace(lower_plot, upper_plot, 800)
    ys = huber_penalty_fun(us).full()

    # plot
    latexify_plot()
    plt.figure(figsize=(10, 10))
    plt.plot(us, ys, label='Huber penalty')

    u_lins = [5, 10.2, 12]
    for i in range(len(u_lins)):
        ys_quad = huber_ggn_quad_fun(us, u_lins[i]).full()
        plt.plot(us, ys_quad, '--', label=f'GGN approximation at {u_lins[i]:.1f}')
        ys_quad = huber_xgn_quad_fun(us, u_lins[i]).full()
        plt.plot(us, ys_quad, ':', label=f'XGN approximation at {u_lins[i]:.1f}')

    plt.grid()
    plt.legend()
    plt.xlim(us[0], us[-1])
    plt.ylim(min(ys), max(ys))
    plt.show()


def plot_one_sided_huber_through_ocp():

    delta = 1e-1
    weight = 1e2
    lower = None
    upper = .0

    # create dummy ocp
    ocp = AcadosOcp()
    ocp.cost.cost_type = "CONVEX_OVER_NONLINEAR"
    ocp.model.cost_y_expr = ca.SX.zeros(0)
    ocp.model.cost_psi_expr = ca.SX.zeros(1)
    ocp.model.cost_r_in_psi_expr = ca.SX.zeros(0)
    x = ca.SX.sym('x', 1)
    ocp.model.x = x
    ocp.formulate_constraint_as_Huber_penalty(x, weight, upper_bound=upper, lower_bound=lower, huber_delta=delta)

    # extract penalty
    penalty_expr = ca.substitute(ocp.model.cost_psi_expr,
                                 ocp.model.cost_r_in_psi_expr,
                                 ocp.model.cost_y_expr - ocp.cost.yref)
    huber_penalty_fun = ca.Function('penalty', [x], [penalty_expr])

    # extract quadratic approximation
    x_ = ca.SX.sym('x_')
    penalty_hess_xgn = ca.substitute(ocp.model.cost_conl_custom_outer_hess,
                                     ocp.model.cost_r_in_psi_expr,
                                     ocp.model.cost_y_expr - ocp.cost.yref)
    penalty_hess_ggn, penalty_grad = ca.hessian(penalty_expr, x)
    # x linearization point, x_ for evaluation
    huber_ggn_quad_fun = ca.Function('penalty', [x_, x], [penalty_expr + penalty_grad*(x_ - x) + 0.5*penalty_hess_ggn*(x_ - x)**2])

    # compute normalized constraint expression
    normalized_constr_expr = x - upper
    J_normalized = ca.jacobian(normalized_constr_expr, x)

    huber_xgn_quad_fun = ca.Function('penalty', [x_, x], [penalty_expr + penalty_grad*(x_ - x) + 0.5*J_normalized.T@penalty_hess_xgn@J_normalized*(x_ - x)**2])

    # create plot data
    upper_plot = 2.0
    lower_plot = -1.
    us = np.linspace(lower_plot, upper_plot, 100)
    ys = huber_penalty_fun(us).full()

    # plot
    latexify_plot()
    plt.figure(figsize=(10, 10))
    plt.plot(us, ys, label='Huber penalty')

    u_lins = np.array([.5, 1.05, 1.5])
    for i in range(len(u_lins)):
        ys_quad = huber_ggn_quad_fun(us, u_lins[i]).full()
        plt.plot(us, ys_quad, '-', alpha = .5, label=f'GGN approximation at {u_lins[i]:.1f}')
        ys_quad = huber_xgn_quad_fun(us, u_lins[i]).full()
        plt.plot(us, ys_quad, '--', linewidth=2, label=f'XGN approximation at {u_lins[i]:.1f}')

    plt.grid()
    plt.legend()
    plt.xlim(us[0], us[-1])
    plt.ylim(min(ys), max(ys))
    plt.show()


if __name__ == "__main__":
    plot_huber_penalty(symmetric=False)
    plot_one_sided_huber_through_ocp()
    # plot_huber_penalty(symmetric=True)
    # plot_huber_through_ocp()

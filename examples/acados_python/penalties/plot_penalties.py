import casadi as ca
import matplotlib.pyplot as plt
import numpy as np

from acados_template import latexify_plot, symmetric_huber_penalty, one_sided_huber_penalty

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


if __name__ == "__main__":
    plot_huber_penalty(symmetric=True)
    plot_huber_penalty(symmetric=False)

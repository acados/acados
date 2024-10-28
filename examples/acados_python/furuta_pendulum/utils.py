import matplotlib.pyplot as plt
from acados_template import latexify_plot
import numpy as np
latexify_plot()

def plot_furuta_pendulum(t_sim, X_sim, U_sim, u_max, plt_show=True):

    nx = 4
    nu = 1
    fig, axes = plt.subplots(nrows=nx+nu, ncols=1, sharex=True)

    labels_x = [r'$\theta_1$', r'$\theta_2$', r'$\dot{\theta}_1$', r'$\dot{\theta}_2$']

    for i, (ax, l) in enumerate(zip(axes, labels_x)):
        ax.plot(t_sim, X_sim[:, i], label=l)
        ax.set_ylabel(l)
        ax.grid()
        ax.set_xlim(t_sim[0], t_sim[-1])

    axes[-1].axhline(u_max, color='k', alpha=0.5, linestyle='dashed')
    axes[-1].axhline(-u_max, color='k', alpha=0.5, linestyle='dashed')
    axes[-1].stairs(U_sim.ravel(), t_sim)
    axes[-1].set_ylabel(r"$\tau_1$")
    axes[-1].grid()
    axes[-1].set_xlabel(r'$t$')

    if plt_show:
        plt.show()


def plot_time_per_solve(times, timeout_max_time: float = 0, heuristic='', plt_show=True, store_figure=False):
    fig, axes = plt.subplots(nrows=1, ncols=1)

    num_solves = times.shape[0]
    if timeout_max_time > 0:
        axes.axhline(timeout_max_time, label='maximum time', color='k', linestyle='dashed', alpha=0.8)
        axes.legend()
        title = f"With timeout using heuristic {heuristic}"
    else:
        title = "No timeout"

    axes.set_title(title)
    axes.bar(np.arange(num_solves), height=times, width=0.5*np.ones((num_solves, )))
    axes.set_xlabel('iteration $k$')
    axes.set_ylabel('total time in ms')
    axes.set_xlim(0, num_solves)
    if plt_show:
        plt.show()
    if store_figure:
        fig.savefig(f"computation_time_{heuristic if timeout_max_time > 0 else 'NO_TIMEOUT'}.png")

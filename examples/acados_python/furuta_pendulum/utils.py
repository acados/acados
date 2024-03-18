import matplotlib.pyplot as plt
from acados_template import latexify_plot

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
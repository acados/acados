from matplotlib import pyplot as plt
from acados_template import latexify_plot

def plot_trajectory(list_trajectories:list, list_labels:list):

    latexify_plot()
    plt.figure(figsize=(4.5, 2.9))
    for i in range(len(list_trajectories)):
        x1 = list_trajectories[i][0,:]
        x2 = list_trajectories[i][1,:]
        plt.plot(x1, x2, label=list_labels[i])
    plt.xlabel("$x_1$", fontsize=14)
    plt.ylabel("$x_2$", fontsize=14)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    # plt.xlim((-0.2, 0.5))
    # plt.ylim((0.0, 0.6))
    plt.show()

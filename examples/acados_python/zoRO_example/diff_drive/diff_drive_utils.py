import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

from mpc_parameters import MPCParam


def get_latex_plot_params():
    params = {'backend': 'ps',
            'text.latex.preamble': r"\usepackage{gensymb} \usepackage{amsmath}",
            'axes.labelsize': 12,
            'axes.titlesize': 12,
            'legend.fontsize': 12,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'text.usetex': True,
            'font.family': 'serif'
    }

    return params


def plot_timings(timing_dict):
    # latexify plot
    params = get_latex_plot_params()
    matplotlib.rcParams.update(params)

    print("timings\t\tmin\tmean\tmax\n--------------------------------")
    for k, v in timing_dict.items():
        print(f"{k:10}\t{np.min(v):.3f}\t {np.mean(v):.3f}\t {np.max(v):.3f}")

    medianprops = dict(linestyle='-', linewidth=2.5, color='darkgreen')
    green_square = dict(markerfacecolor='palegreen', marker='D')
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    ax.boxplot(timing_dict.values(), vert=False,
        flierprops=green_square, \
        medianprops=medianprops, showmeans=False)
    ax.set_yticklabels(timing_dict.keys())
    plt.grid()
    plt.xlabel('CPU time [ms]')
    plt.tight_layout()
    plt.savefig(os.path.join("figures", "timings_diff_drive.pdf"),
        bbox_inches='tight', transparent=True, pad_inches=0.05)


def compute_min_dis(cfg:MPCParam, s:np.ndarray)->float:
    min_dist = np.inf
    for idx_obs in range(cfg.num_obs):
        min_dist = np.min([min_dist, \
            np.linalg.norm(s[:2] - cfg.obs_pos[idx_obs,:])-cfg.obs_radius[idx_obs]])
    return min_dist
from tokenize import Double
import numpy as np
import matplotlib.pyplot as plt

from mpc_parameters import MPCParam


def plot_timings(timing_dict):
    print("timings\t\tmin\tmean\tmax\n--------------------------------")
    for k, v in timing_dict.items():
        print(f"{k:10}\t{np.min(v):.3f}\t {np.mean(v):.3f}\t {np.max(v):.3f}")

    medianprops = dict(linestyle='-', linewidth=2.5, color='darkgreen')
    green_square = dict(markerfacecolor='palegreen', marker='D')
    plt.rcParams.update({'font.size': 16})
    _, ax = plt.subplots()
    ax.boxplot(timing_dict.values(), vert=False,
        # flierprops=green_square, \
        medianprops=medianprops, showmeans=False)
    ax.set_yticklabels(timing_dict.keys())
    plt.grid()
    plt.xlabel('CPU time [ms]')
    plt.tight_layout()
    plt.savefig("figures/timings_diff_drive.png",
        bbox_inches='tight', transparent=True, pad_inches=0.05)
    plt.show()


def compute_min_dis(cfg:MPCParam, s:np.ndarray)->Double:
    min_dist = np.inf
    for idx_obs in range(cfg.num_obs):
        min_dist = np.min([min_dist, \
            np.linalg.norm(s[:2] - cfg.obs_pos[idx_obs,:])-cfg.obs_radius[idx_obs]])
    return min_dist
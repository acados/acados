from tokenize import Double
import numpy as np
import matplotlib.pyplot as plt

from mpc_parameters import MPCParam


def samplesFromEllipsoid(N, w, Z)->np.ndarray:
    """
    draws samples from ellipsoid with center w and variability matrix Z
    """

    nw = w.shape[0]                  # dimension
    lam, v = np.linalg.eig(Z)

    # sample in hypersphere
    r = np.random.rand()**(1/nw)     # radial position of sample
    x = np.random.randn(N, nw)
    y = np.zeros((N, nw))
    for idx in range(N):
        x[idx,:] = x[idx,:] / np.linalg.norm(x[idx,:])
        x[idx,:] *= r
        # project to ellipsoid
        y[idx,:] = v @ (np.sqrt(lam) * x[idx,:]) + w

    return y


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
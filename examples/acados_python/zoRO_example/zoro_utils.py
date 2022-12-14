import numpy as np


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
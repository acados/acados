from dataclasses import dataclass, field
import numpy as np


@dataclass
class ZoroDescription:
    """
    Zero-Order Robust Optimization scheme.

    For advanced users.
    """
    backoff_scaling_gamma: float = 1.0
    fdbk_K_mat: np.ndarray = None
    unc_jac_G_mat: np.ndarray = None    # default: an identity matrix
    P0_mat: np.ndarray = None
    W_mat: np.ndarray = None
    idx_lbx_t: list = field(default_factory=list)
    idx_ubx_t: list = field(default_factory=list)
    idx_lbx_e_t: list = field(default_factory=list)
    idx_ubx_e_t: list = field(default_factory=list)
    idx_lbu_t: list = field(default_factory=list)
    idx_ubu_t: list = field(default_factory=list)
    idx_lg_t: list = field(default_factory=list)
    idx_ug_t: list = field(default_factory=list)
    idx_lg_e_t: list = field(default_factory=list)
    idx_ug_e_t: list = field(default_factory=list)
    idx_lh_t: list = field(default_factory=list)
    idx_uh_t: list = field(default_factory=list)
    idx_lh_e_t: list = field(default_factory=list)
    idx_uh_e_t: list = field(default_factory=list)

def process_zoro_description(zoro_description: ZoroDescription):
    zoro_description.nw, _ = zoro_description.W_mat.shape
    if zoro_description.unc_jac_G_mat is None:
        zoro_description.unc_jac_G_mat = np.eye(zoro_description.nw)
    zoro_description.nlbx_t = len(zoro_description.idx_lbx_t)
    zoro_description.nubx_t = len(zoro_description.idx_ubx_t)
    zoro_description.nlbx_e_t = len(zoro_description.idx_lbx_e_t)
    zoro_description.nubx_e_t = len(zoro_description.idx_ubx_e_t)
    zoro_description.nlbu_t = len(zoro_description.idx_lbu_t)
    zoro_description.nubu_t = len(zoro_description.idx_ubu_t)
    zoro_description.nlg_t = len(zoro_description.idx_lg_t)
    zoro_description.nug_t = len(zoro_description.idx_ug_t)
    zoro_description.nlg_e_t = len(zoro_description.idx_lg_e_t)
    zoro_description.nug_e_t = len(zoro_description.idx_ug_e_t)
    zoro_description.nlh_t = len(zoro_description.idx_lh_t)
    zoro_description.nuh_t = len(zoro_description.idx_uh_t)
    zoro_description.nlh_e_t = len(zoro_description.idx_lh_e_t)
    zoro_description.nuh_e_t = len(zoro_description.idx_uh_e_t)
    return zoro_description.__dict__

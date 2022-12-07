from dataclasses import dataclass, field
import numpy as np


@dataclass
class ZoroDescription:
    fdbk_K_mat: np.ndarray = None   # TODO: default: a full-zero matrix
    unc_jac_G_mat: np.ndarray = None    # default: an identity matrix
    P0_mat: np.ndarray = None       # TODO: default: a full-zero matrix
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

def process_zoro_stuff(zoro_stuff: ZoroDescription):
    zoro_stuff.nw, _ = zoro_stuff.W_mat.shape
    if zoro_stuff.unc_jac_G_mat is None:
        zoro_stuff.unc_jac_G_mat = np.eye(zoro_stuff.nw)
    zoro_stuff.nlbx_t = len(zoro_stuff.idx_lbx_t)
    zoro_stuff.nubx_t = len(zoro_stuff.idx_ubx_t)
    zoro_stuff.nlbx_e_t = len(zoro_stuff.idx_lbx_e_t)
    zoro_stuff.nubx_e_t = len(zoro_stuff.idx_ubx_e_t)
    zoro_stuff.nlbu_t = len(zoro_stuff.idx_lbu_t)
    zoro_stuff.nubu_t = len(zoro_stuff.idx_ubu_t)
    zoro_stuff.nlg_t = len(zoro_stuff.idx_lg_t)
    zoro_stuff.nug_t = len(zoro_stuff.idx_ug_t)
    zoro_stuff.nlg_e_t = len(zoro_stuff.idx_lg_e_t)
    zoro_stuff.nug_e_t = len(zoro_stuff.idx_ug_e_t)
    zoro_stuff.nlh_t = len(zoro_stuff.idx_lh_t)
    zoro_stuff.nuh_t = len(zoro_stuff.idx_uh_t)
    zoro_stuff.nlh_e_t = len(zoro_stuff.idx_lh_e_t)
    zoro_stuff.nuh_e_t = len(zoro_stuff.idx_uh_e_t)
    return zoro_stuff.__dict__
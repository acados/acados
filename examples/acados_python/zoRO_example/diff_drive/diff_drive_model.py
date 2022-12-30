from enum import Enum
import numpy as np

from acados_template import AcadosModel
from casadi import SX, vertcat, cos, sin, sqrt, sumsqr

from mpc_parameters import MPCParam


class RobotState(Enum):
    POSX = 0
    POSY = 1
    THETA = 2
    VEL = 3
    OMEGA = 4


def dxdt(x0, u):
    dx = vertcat(x0[3]*cos(x0[2]),
                 x0[3]*sin(x0[2]),
                 x0[4],
                 u[0],
                 u[1])
    return dx


def export_diff_drive_model(as_cfg:MPCParam):

    model_name = 'diff_drive_model'

    nx = as_cfg.nx
    nu = as_cfg.nu
    num_obs = as_cfg.num_obs

    x = SX.sym('x', nx, 1)
    u = SX.sym('u', nu, 1)
    xdot = SX.sym('xdot', nx, 1)

    # dynamics
    f_expl = vertcat(x[3]*cos(x[2]),
                     x[3]*sin(x[2]),
                     x[4],
                     u[0],
                     u[1])

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_expl_expr = f_expl
    model.f_impl_expr = f_impl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = model_name

    con_expr = SX.sym('con_expr', num_obs, 1)
    obs_param = SX.sym('obs', 2*num_obs, 1)
    # [x0, y0, x1, y1, ...]
    p_obs = [obs_param[2*k:2*(k+1)] for k in range(0, num_obs)]

    for i in range(0, num_obs):
        dist = x[0:2] - p_obs[i]
        con_expr[i] = sqrt(sumsqr(dist))
    model.p = obs_param
    model.con_h_expr = con_expr
    model.con_h_expr_e = con_expr

    return model
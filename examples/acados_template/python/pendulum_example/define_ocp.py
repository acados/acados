from acados_template import *
import acados_template as at
from export_ode_model import *
import numpy as np
import scipy.linalg
from ctypes import *

def define_ocp(model, acados_path='/usr/lib'):
    # create render arguments
    ocp = acados_ocp_nlp()

    # set model_name
    ocp.model_name = model.name

    Tf = 2.0
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx
    N = 50

    # set ocp_nlp_dimensions
    nlp_dims     = ocp.dims
    nlp_dims.nx  = nx
    nlp_dims.ny  = ny
    nlp_dims.ny_e = ny_e
    nlp_dims.nbx = 0
    nlp_dims.nbu = nu
    nlp_dims.nu  = model.u.size()[0]
    nlp_dims.N   = N

    # set weighting matrices
    nlp_cost = ocp.cost
    Q = np.eye(4)
    Q[0,0] = 1e0
    Q[1,1] = 1e2
    Q[2,2] = 1e-3
    Q[3,3] = 1e-2

    R = np.eye(1)
    R[0,0] = 1e0

    nlp_cost.W = scipy.linalg.block_diag(Q, R)

    Vx = np.zeros((ny, nx))
    Vx[0,0] = 1.0
    Vx[1,1] = 1.0
    Vx[2,2] = 1.0
    Vx[3,3] = 1.0

    nlp_cost.Vx = Vx

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    nlp_cost.Vu = Vu

    nlp_cost.W_e = Q

    Vx_e = np.zeros((ny_e, nx))
    Vx_e[0,0] = 1.0
    Vx_e[1,1] = 1.0
    Vx_e[2,2] = 1.0
    Vx_e[3,3] = 1.0

    nlp_cost.Vx_e = Vx_e

    nlp_cost.yref  = np.zeros((ny, ))
    nlp_cost.yref_e = np.zeros((ny_e, ))

    # setting bounds
    Fmax = 2.0
    nlp_con = ocp.constraints
    nlp_con.lbu = np.array([-Fmax])
    nlp_con.ubu = np.array([+Fmax])
    nlp_con.x0 = np.array([0.0, 3.14, 0.0, 0.0])
    # nlp_con.x0 = np.array([0.0, 0.5, 0.0, 0.0])
    nlp_con.idxbu = np.array([0])

    # set constants
    # ocp.constants['PI'] = 3.1415926535897932

    # set QP solver
    # ocp.solver_config.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_config.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_config.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_config.integrator_type = 'ERK'

    # set prediction horizon
    ocp.solver_config.tf = Tf
    ocp.solver_config.nlp_solver_type = 'SQP'
    # ocp.solver_config.nlp_solver_type = 'SQP_RTI'

    # set header path
    ocp.acados_include_path  = f'{acados_path}/include'
    ocp.acados_lib_path      = f'{acados_path}/lib'
    return ocp

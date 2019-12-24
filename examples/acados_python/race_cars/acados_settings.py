#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

# author: Daniel Kloeser

from acados_template import *
from bycicle_model import *
import scipy.linalg
import numpy as np


def acados_settings(Tf,N,track_file):
    # create render arguments
    ocp = acados_ocp_nlp()

    # export model
    model,constraint = bycicle_model(track_file)

    # define acados ODE
    model_ac = acados_dae()
    model_ac.f_impl_expr = model.f_impl_expr
    model_ac.f_expl_expr = model.f_expl_expr
    model_ac.x = model.x
    model_ac.xdot = model.xdot
    model_ac.u = model.u
    model_ac.z = model.z
    model_ac.p = model.p
    model_ac.name = model.name
    ocp.model=model_ac

    # define constraint
    constraint_ac=acados_constraint()
    constraint_ac.con_h_expr=constraint.expr
    constraint_ac.x=constraint.x
    constraint_ac.u=constraint.u
    constraint_ac.p=constraint.p
    constraint_ac.nh=constraint.nc
    constraint_ac.z=vertcat([])
    constraint_ac.name="con_bycicle"

    # set ocp_nlp_dimensions
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx
    nlp_dims     = ocp.dims
    nlp_dims.nx  = nx
    nlp_dims.np = 0
    nlp_dims.ny  = ny
    nlp_dims.ny_e = ny_e
    nlp_dims.nbx = 1
    nlp_dims.nsbx = 0
    nlp_dims.nbu = nu
    nlp_dims.nu  = nu
    nlp_dims.N   = N
    nlp_dims.nsh=2
    nlp_dims.nh=constraint.nc
    nlp_dims.ns=2

    # set weighting matrices
    Q=np.zeros((nx,nx))
    Q[0, 0] = 1e-1
    Q[1, 1] = 1e-8
    Q[2, 2] = 1e-8
    Q[3, 3] = 1e-8
    Q[4, 4] = 1e-3
    Q[5, 5] = 5e-3

    R=np.eye(nu)
    R[0, 0] = 1e-3
    R[1, 1] = 5e-3

    Qe=np.zeros((nx,nx))
    Qe[0, 0] = 5e0
    Qe[1, 1] = 1e1
    Qe[2, 2] = 1e-8
    Qe[3, 3] = 1e-8
    Qe[4, 4] = 5e-3
    Qe[5, 5] = 2e-3

    nlp_cost = ocp.cost
    nlp_cost.cost_type = 'LINEAR_LS'
    nlp_cost.cost_type_e = 'LINEAR_LS'
    nlp_cost.W = scipy.linalg.block_diag(Q, R)
    nlp_cost.W_e = Qe
    
    Vx = np.zeros((ny, nx))
    Vx[0,0] = 1.0
    Vx[1,1] = 1.0
    Vx[2,2] = 1.0
    Vx[3,3] = 1.0
    Vx[4,4] = 1.0
    Vx[5,5] = 1.0
    nlp_cost.Vx = Vx

    Vu = np.zeros((ny, nu))
    Vu[6,0] = 1.0
    Vu[7,1] = 1.0
    nlp_cost.Vu = Vu

    Vx_e = np.zeros((ny_e, nx))
    Vx_e[0,0] = 1.0
    Vx_e[1,1] = 1.0
    Vx_e[2,2] = 1.0
    Vx_e[3,3] = 1.0
    Vx_e[4,4] = 1.0
    Vx_e[5,5] = 1.0
    nlp_cost.Vx_e = Vx_e

    nlp_cost.zl = 100 * np.ones((nlp_dims.ns,))
    nlp_cost.Zl = 0 * np.ones((nlp_dims.ns,))
    nlp_cost.zu = 100 * np.ones((nlp_dims.ns,))
    nlp_cost.Zu = 0 * np.ones((nlp_dims.ns,))

    # set intial references
    nlp_cost.yref  = np.array([1,0,0,0,0,0,0,0])
    nlp_cost.yref_e = np.array([0,0,0,0,0,0])

    # setting bounds
    ocp.con_h=constraint_ac
    nlp_con = ocp.constraints
    nlp_con.lbx=np.array([-12])
    nlp_con.ubx=np.array([12])
    nlp_con.idxbx=np.array([1])
    nlp_con.lbu = np.array([model.dthrottle_min,model.ddelta_min])
    nlp_con.ubu = np.array([model.dthrottle_max,  model.ddelta_max])
    nlp_con.idxbu = np.array([0,1])
    #nlp_con.lsbx=np.zero s([1])
    #nlp_con.usbx=np.zeros([1])
    #nlp_con.idxsbx=np.array([1])
    nlp_con.lh=np.array([constraint.along_min,constraint.alat_min,model.n_min,model.throttle_min,model.delta_min])
    nlp_con.uh=np.array([constraint.along_max,constraint.alat_max,model.n_max,model.throttle_max,model.delta_max])
    nlp_con.lsh=np.zeros(nlp_dims.nsh)
    nlp_con.ush=np.zeros(nlp_dims.nsh)
    nlp_con.idxsh=np.array([0,2])

    # set intial condition
    nlp_con.x0 = model.x0

    # set QP solver and integration
    ocp.solver_options.tf = Tf
    # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.sim_method_num_stages=4
    ocp.solver_options.sim_method_num_steps=3

    # ocp.solver_options.qp_solver_tol_stat = 1e-2
    # ocp.solver_options.qp_solver_tol_eq = 1e-2
    # ocp.solver_options.qp_solver_tol_ineq = 1e-2
    # ocp.solver_options.qp_solver_tol_comp = 1e-2

    # create solver
    acados_solver = generate_solver(ocp, json_file='acados_ocp.json')

    return constraint,model,acados_solver

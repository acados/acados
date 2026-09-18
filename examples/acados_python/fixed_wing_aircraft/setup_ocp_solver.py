#
# Copyright (c) The acados authors.
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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
import casadi as ca
import numpy as np
import scipy.linalg as la

from model import export_aircraft_model


def setup_ocp_solver(x0, N_horizon, Tf, RTI=True, soft_constraints=True):
    ocp = AcadosOcp()

    # set model
    model = export_aircraft_model()
    ocp.model = model

    # set parameters (wind)
    ocp.parameter_values = np.zeros(3)  # default: no wind

    nx = model.x.rows()
    nu = model.u.rows()

    ### cost ###
    ocp.cost.cost_type   = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    x_track = model.x

    ocp.model.cost_y_expr   = ca.vertcat(x_track, model.u)
    ocp.model.cost_y_expr_e = x_track

    ny   = x_track.rows() + nu
    ny_e = x_track.rows()

    Q_mat = np.diag([.1]*3 + [.1]*3 + [1]*9 + [1]*3 + [0.1]*2 + [1]*2)
    R_mat = np.diag([.1, .1, 50, 50])

    ocp.cost.W   = la.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    ocp.cost.yref   = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    ### constraints ###
    # control constraints
    ocp.constraints.lbu = np.array([-5, -5, np.deg2rad(-19), np.deg2rad(-19)])
    ocp.constraints.ubu = np.array([ 5,  5, np.deg2rad( 19), np.deg2rad( 19)])
    ocp.constraints.idxbu = np.array([0, 1, 2, 3])

    # state constraints
    ocp.constraints.lbx = np.array([-1, -1, -1,  0,  0, np.deg2rad(-32), np.deg2rad(-25)])
    ocp.constraints.ubx = np.array([ 1,  1,  1, 20, 20, np.deg2rad( 32), np.deg2rad( 25)])
    ocp.constraints.idxbx = np.array([15, 16, 17, 18, 19, 20, 21])
    ocp.constraints.x0 = x0

    # nonlinear constraints (body velocity, angle of attack, sideslip angle)
    ocp.dims.nh = model.con_h_expr.rows()
    ocp.constraints.lh = np.array([12, -5, -5, np.deg2rad(-6), np.deg2rad(-10)])
    ocp.constraints.uh = np.array([20,  5,  5, np.deg2rad(12), np.deg2rad( 10)])

    # prediction horizon
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = Tf

    # soft constraints
    if soft_constraints:
        ocp.constraints.idxsh = np.arange(5)

        # linear and quadratic slack cost
        ocp.cost.zl = 1e2*np.ones(5)
        ocp.cost.zu = 1e2*np.ones(5)
        ocp.cost.Zl = 1e1*np.ones(5)
        ocp.cost.Zu = 1e1*np.ones(5)

        ocp.constraints.idxsh_e = ocp.constraints.idxsh
        ocp.model.con_h_expr_e = ocp.model.con_h_expr
        ocp.constraints.lh_e = ocp.constraints.lh
        ocp.constraints.uh_e = ocp.constraints.uh
        ocp.cost.zl_e = ocp.cost.zl
        ocp.cost.zu_e = ocp.cost.zu
        ocp.cost.Zl_e = ocp.cost.Zl
        ocp.cost.Zu_e = ocp.cost.Zu

    if RTI:
        ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    else:
        ocp.solver_options.nlp_solver_type = 'SQP'

    ### ocp solver options ###
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.collocation_type = 'GAUSS_RADAU_IIA'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.levenberg_marquardt = 1e-4
    ocp.solver_options.sim_method_newton_iter = 10

    ocp.code_export_directory = 'c_generated_code/ocp'
    ocp_solver = AcadosOcpSolver(ocp, json_file='data/ocp_solver.json')

    return ocp, ocp_solver


def setup_integrator(ocp):
    sim = AcadosSim.from_ocp(ocp)

    sim.parameter_values = np.zeros(3)  # default: no wind

    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.collocation_type = 'GAUSS_RADAU_IIA'
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 4
    sim.solver_options.newton_iter = 3

    # needed for the EKF
    sim.solver_options.sens_forw = True
    sim.solver_options.sens_forw_p = True

    sim.code_export_directory = 'c_generated_code/sim'

    return AcadosSimSolver(sim, json_file='data/simulator.json')

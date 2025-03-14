# -*- coding: future_fstrings -*-
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

import numpy as np
import casadi as ca

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosMultiphaseOcp, AcadosModel, ACADOS_INFTY


def export_double_integrator_model(dim_q, dt) -> AcadosModel:
    model_name = 'double_integrator_disc_dyn'

    # set up states & controls
    q = ca.SX.sym('q', dim_q)
    v = ca.SX.sym('v', dim_q)

    x = ca.vertcat(q, v)

    u = ca.SX.sym('u', dim_q)

    vnext = v + u*dt
    qnext = q + vnext*dt

    model = AcadosModel()

    model.disc_dyn_expr = ca.vertcat(qnext, vnext)
    model.x = x
    model.u = u
    model.name = model_name

    return model

def main():
    # Horizon definition
    N = 20 # Number of time intervals
    Tf = 1.0 # Duration
    dt = Tf/N # Time step

    # Dimension of double integrator configuration and velocity
    nq = 1
    nv = 1

    # State and control dimensions
    nx = nq + nv
    nu = nv

    # Initial state
    x0 = np.concatenate((np.ones(nq), 0.25*np.ones(nv)))

    # Cost matrices
    Q_mat = np.diag(np.concatenate((np.ones(nq), 10*np.ones(nv))))
    R_mat = np.diag(10*np.ones(nu))

    # Number of intervals in each phase
    N_list = [10, 10]
    # Multiphase ocp
    multiphase_ocp = AcadosMultiphaseOcp(N_list=N_list)

    #### PHASE 0: velocity can take any value ####
    # Dynamics
    phase_idx = 0
    acados_model = export_double_integrator_model(nq, dt)

    ocp = AcadosOcp()
    ocp.model = acados_model
    ocp.solver_options.N_horizon = N_list[phase_idx]

    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    # Quadratic state and control cost
    ocp.model.cost_expr_ext_cost = acados_model.x.T @ Q_mat @ acados_model.x + acados_model.u.T @ R_mat @ acados_model.u
    ocp.model.cost_expr_ext_cost_e = acados_model.x.T @ Q_mat @ acados_model.x

    # Control limits
    Fmax = 1
    ocp.constraints.lbu = -Fmax*np.ones(nu)
    ocp.constraints.ubu = Fmax*np.ones(nu)
    ocp.constraints.idxbu = np.arange(nu)

    # Initial state constraint
    ocp.constraints.x0 = x0

    multiphase_ocp.set_phase(ocp, phase_idx)

    #### PHASE 1: velocity must be nonpositive ####
    phase_idx += 1

    # Dynamics
    acados_model = export_double_integrator_model(nq, dt)

    # Nonlinear constraint: velocity must be nonpositive (yes, this is actually
    # a linear constraint, but we're enforcing it using the acados nonlinear constraint interface)
    h_expr = acados_model.x[nq:]
    h_lb = -1e-1 * ACADOS_INFTY*np.ones(nv)  # One-sided constraint with qpOASES not supported yet, otherwise we'd use -ACADOS_INFTY
    h_ub = np.zeros(nv)

    acados_model.con_h_expr = h_expr

    ocp = AcadosOcp()
    ocp.model = acados_model
    ocp.solver_options.N_horizon = N_list[phase_idx]

    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    # Quadratic state and control cost
    ocp.model.cost_expr_ext_cost = acados_model.x.T @ Q_mat @ acados_model.x + acados_model.u.T @ R_mat @ acados_model.u
    ocp.model.cost_expr_ext_cost_e = acados_model.x.T @ Q_mat @ acados_model.x

    # Control limits
    Fmax = 1
    ocp.constraints.lbu = -Fmax*np.ones(nu)
    ocp.constraints.ubu = Fmax*np.ones(nu)
    ocp.constraints.idxbu = np.arange(nu)

    # Set nonlinear constraint sizes and bounds
    ocp.dims.nh = h_expr.shape[0]
    ocp.constraints.lh = h_lb
    ocp.constraints.uh = h_ub

    multiphase_ocp.set_phase(ocp, phase_idx)

    # Set options
    multiphase_ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    multiphase_ocp.solver_options.qp_solver_cond_N = N
    multiphase_ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    multiphase_ocp.solver_options.nlp_solver_type = 'SQP'
    multiphase_ocp.solver_options.tf = Tf
    multiphase_ocp.solver_options.nlp_solver_tol_eq = 1e-4
    multiphase_ocp.solver_options.nlp_solver_tol_ineq = 1e-4

    multiphase_ocp.mocp_opts.integrator_type = ['DISCRETE', 'DISCRETE']

    ocp_solver = AcadosOcpSolver(multiphase_ocp, json_file = 'acados_ocp.json')

    status = ocp_solver.solve()
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
    assert status == 0, f'acados returned status {status}'

if __name__ == '__main__':
    main()

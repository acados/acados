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

import casadi as ca
import numpy as np
from acados_template import AcadosOcp, AcadosOcpSolver, ocp_get_default_cmake_builder
import scipy.linalg

def linear_mass_spring_model(casadi_dynamics, casadi_cost):

    # dims
    num_mass = 4

    nx = 2*num_mass
    nu = num_mass-1

    # symbolic variables
    x = ca.SX.sym('x', nx, 1) # states
    u = ca.SX.sym('u', nu, 1) # controls

    # dynamics
    # continuous time
    Ac = np.zeros((nx, nx))
    for ii in range(num_mass):
        Ac[ii,num_mass+ii] = 1.0
        Ac[num_mass+ii,ii] = -2.0

    for ii in range(num_mass-1):
        Ac[num_mass+ii,ii+1] = 1.0
        Ac[num_mass+ii+1,ii] = 1.0


    Bc = np.zeros((nx, nu))
    for ii in range(nu):
        Bc[num_mass+ii, ii] = 1.0

    c_const = np.zeros(nx)

    # discrete time
    Ts = 0.5 # sampling time
    M = scipy.linalg.expm(np.vstack(( np.hstack((Ts*Ac, Ts*Bc)), np.zeros((nu, int(2*nx/2+nu))))))
    A = M[:nx,:nx]
    B = M[:nx,nx:]

    expr_f_expl = Ac @ x + Bc @ u + c_const
    discrete_dynamics_expr = A @ x + B @ u

    # constraints
    expr_h = ca.vertcat(u, x)
    expr_h_e = x

    # external cost
    yr_u = np.zeros(nu)
    yr_x = np.zeros(nx)
    dWu = 2*np.ones(nu)
    dWx = np.ones(nx)

    ymyr = ca.vertcat(u, x) - ca.vertcat(yr_u, yr_x)
    ymyr_e = x - yr_x

    expr_ext_cost = 0.5 * ymyr.T @ (np.concatenate((dWu, dWx)) *  ymyr)
    expr_ext_cost_e = 0.5 * ymyr_e.T @ (dWx * ymyr_e)

    x0 = np.zeros(nx)
    x0[0] = 2.5
    x0[1] = 2.5

    lh = - np.concatenate(( 0.5 * np.ones(nu), 4.0 * np.ones(nx)))
    uh = + np.concatenate(( 0.5 * np.ones(nu), 4.0 * np.ones(nx)))
    lh_e = -4.0 * np.ones(nx)
    uh_e = 4.0 * np.ones(nx)

    # acados ocp model
    ocp = AcadosOcp()
    ocp.model.name = 'lin_mass'

    # symbolics
    ocp.model.x = x
    ocp.model.u = u

    # cost
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    if not casadi_dynamics:
        # Generic dynamics
        ocp.model.dyn_ext_fun_type = 'generic'
        ocp.model.dyn_generic_source = 'generic_disc_dyn.c'
        ocp.model.dyn_disc_fun = 'disc_dyn_fun'
        ocp.model.dyn_disc_fun_jac = 'disc_dyn_fun_jac'
        ocp.model.dyn_disc_fun_jac_hess = 'disc_dyn_fun_jac_hess' # only needed for exact hessi
    else:
        # dynamics expression
        ocp.model.disc_dyn_expr = discrete_dynamics_expr

    if not casadi_cost:
        # Generic stage cost
        ocp.model.cost_ext_fun_type = 'generic'
        ocp.model.cost_source_ext_cost = 'generic_ext_cost.c'
        ocp.model.cost_function_ext_cost = 'ext_cost'
        # Generic terminal cost
        ocp.model.cost_ext_fun_type_e = 'generic'
        ocp.model.cost_source_ext_cost_e = 'generic_ext_cost.c'
        ocp.model.cost_function_ext_cost_e = 'ext_costN'
    else:
        # cost expression
        ocp.model.cost_expr_ext_cost = expr_ext_cost
        ocp.model.cost_expr_ext_cost_e = expr_ext_cost_e

    # constraints
    ocp.constraints.x0 = x0

    ocp.model.con_h_expr = expr_h
    ocp.constraints.lh = lh
    ocp.constraints.uh = uh

    ocp.model.con_h_expr_e = expr_h_e
    ocp.constraints.lh_e = lh_e
    ocp.constraints.uh_e = uh_e

    # acados ocp opts
    N = 20
    shooting_nodes = np.linspace(0,10,N+1)
    T = 10.0 # horizon length time

    ocp.solver_options.tf = T
    ocp.solver_options.N_horizon = N
    ocp.solver_options.shooting_nodes = shooting_nodes
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.nlp_solver_ext_qp_res = 0
    ocp.solver_options.nlp_solver_max_iter = 100
    ocp.solver_options.tol = 1e-10
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.qp_solver_cond_N = 5
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.print_level = 2

    return ocp



def main(casadi_dynamics, casadi_cost, use_cmake = False):

    ocp = linear_mass_spring_model(casadi_dynamics, casadi_cost)

    # create ocp solver
    cmake_builder = ocp_get_default_cmake_builder() if use_cmake else None
    ocp_solver = AcadosOcpSolver(ocp, cmake_builder=cmake_builder)

    # solve
    status = ocp_solver.solve()

    # get solution
    u0 = ocp_solver.get(0, 'u')
    x0 = ocp_solver.get(0, 'x')

    # get info
    sqp_iter = ocp_solver.get_stats('sqp_iter')
    time_tot = ocp_solver.get_stats('time_tot')
    time_lin = ocp_solver.get_stats('time_lin')
    time_reg = ocp_solver.get_stats('time_reg')
    time_qp = ocp_solver.get_stats('time_qp')

    print(f'\nstatus = {status}, sqp_iter = {sqp_iter}')

    # print statistics
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception('ocp_nlp solver returned status nonzero')
    elif sqp_iter > 2:
        raise Exception('ocp can be solved in 2 iterations!')
    else:
        print(f'test_ocp_linear_mass_spring: success')

    ocp_solver.store_iterate(filename=f'{ocp.model.name}_{casadi_dynamics}.json', overwrite=True)


if __name__ == '__main__':

    casadi_dynamics = False # False = generic, True = casadi
    casadi_cost = True # False = generic, True = casadi # NOTE: generic cost is not implemented in python interface.

    main(casadi_dynamics, casadi_cost)

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
from acados_template import AcadosOcp, AcadosOcpSolver, casadi_length
import scipy.linalg

def linear_mass_spring_model():

    # dims
    num_mass = 4

    nx = 2*num_mass
    nu = num_mass-1

    # symbolic variables
    sym_x = ca.SX.sym('x', nx, 1) # states
    sym_u = ca.SX.sym('u', nu, 1) # controls
    sym_xdot = ca.SX.sym('xdot', casadi_length(sym_x)) #state derivatives

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

    expr_f_expl = Ac@sym_x + Bc@sym_u + c_const
    expr_f_impl = expr_f_expl - sym_xdot
    expr_phi = A@sym_x + B@sym_u

    # constraints
    expr_h = ca.vertcat(sym_u, sym_x)
    expr_h_e = sym_x

    # nonlnear least squares
    expr_y = ca.vertcat(sym_u, sym_x)
    expr_y_e = sym_x

    # external cost
    yr_u = np.zeros(nu)
    yr_x = np.zeros(nx)
    dWu = 2*np.ones(nu)
    dWx = np.ones(nx)

    ymyr = ca.vertcat(sym_u, sym_x) - ca.vertcat(yr_u, yr_x)
    ymyr_e = sym_x - yr_x

    expr_ext_cost = 0.5 * (ymyr.T @ ( np.concatenate((dWu, dWx)) *  ymyr))
    expr_ext_cost_e = 0.5 * ymyr_e.T @ (dWx * ymyr_e)

    # populate structure
    model = {}
    model['nx'] = nx
    model['nu'] = nu
    model['sym_x'] = sym_x
    model['sym_xdot'] = sym_xdot
    model['sym_u'] = sym_u
    model['expr_f_expl'] = expr_f_expl
    model['expr_f_impl'] = expr_f_impl
    model['expr_phi'] = expr_phi
    model['expr_h'] = expr_h
    model['expr_h_e'] = expr_h_e
    model['expr_y'] = expr_y
    model['expr_y_e'] = expr_y_e
    model['expr_ext_cost'] = expr_ext_cost
    model['expr_ext_cost_e'] = expr_ext_cost_e

    return model



def main():
    N = 20
    tol = 1e-10
    shooting_nodes = np.linspace(0,10,N+1)

    model_name = 'lin_mass'
    nlp_solver = 'SQP'
    # nlp_solver_exact_hessian = 'true'
    regularize_method = 'CONVEXIFY'
    nlp_solver_max_iter = 100
    qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    qp_solver_cond_N = 5
    cost_type = 'EXTERNAL'

    model = linear_mass_spring_model()

    T = 10.0 # horizon length time
    nx = model['nx']
    nu = model['nu']

    x0 = np.zeros(nx)
    x0[0] = 2.5
    x0[1]=2.5

    lh = - np.concatenate(( 0.5 * np.ones(nu), 4.0 * np.ones(nx)))
    uh = + np.concatenate(( 0.5 * np.ones(nu), 4.0 * np.ones(nx)))
    lh_e = -4.0 * np.ones(nx)
    uh_e = 4.0 * np.ones(nx)

    # acados ocp model
    casadi_dynamics = 0 # 0=generic, 1=casadi
    casadi_cost = 1 # 0=generic, 1=casadi # NOTE: generic cost is not implemented in python interface.

    ocp = AcadosOcp()
    ocp.model.name = model_name
    ocp.solver_options.tf = T

    # symbolics
    ocp.model.x = model['sym_x']
    ocp.model.u = model['sym_u']
    ocp.model.xdot = model['sym_xdot']

    # cost
    ocp.cost.cost_type = cost_type
    ocp.cost.cost_type_e = cost_type

    if (casadi_dynamics == 0):
        # Generic dynamics
        ocp.model.dyn_ext_fun_type = 'generic'
        ocp.model.dyn_generic_source = 'generic_disc_dyn.c'
        ocp.model.dyn_disc_fun = 'disc_dyn_fun'
        ocp.model.dyn_disc_fun_jac = 'disc_dyn_fun_jac'
        ocp.model.dyn_disc_fun_jac_hess = 'disc_dyn_fun_jac_hess' # only needed for exact hessi
    else:
        # dynamics expression
        ocp.model.disc_dyn_expr = model['expr_phi']

    if (casadi_cost == 0):
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
        ocp.model.cost_expr_ext_cost = model['expr_ext_cost']
        ocp.model.cost_expr_ext_cost_e = model['expr_ext_cost_e']

    # constraints
    ocp.constraints.x0 = x0
    ocp.model.con_h_expr = model['expr_h']
    ocp.constraints.lh = lh
    ocp.constraints.uh = uh
    ocp.model.con_h_expr_e = model['expr_h_e']
    ocp.constraints.lh_e = lh_e
    ocp.constraints.uh_e = uh_e

    # acados ocp opts
    ocp.dims.N = N
    ocp.solver_options.shooting_nodes = shooting_nodes
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.regularize_method = regularize_method
    ocp.solver_options.nlp_solver_ext_qp_res = 0
    ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter
    ocp.solver_options.tol = tol
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.nlp_solver_type = nlp_solver
    ocp.solver_options.print_level = 2

    # create ocp solver
    ocp_solver = AcadosOcpSolver(ocp)

    # initial state
    ocp_solver.set(0, 'lbx', x0)
    ocp_solver.set(0, 'ubx', x0)

    # initialize
    for i in range(N):
        ocp_solver.set(i, 'x', np.zeros(nx))
        ocp_solver.set(i, 'u', np.zeros(nu))
    ocp_solver.set(N, 'x', np.zeros(nx))

    # solve
    # tic
    status = ocp_solver.solve()
    # time_ext = toc

    # get solution
    utraj = ocp_solver.get(0, 'u')
    xtraj = ocp_solver.get(0, 'x')

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

    ocp_solver.store_iterate(filename=f'{model_name}_{casadi_dynamics}.json', overwrite=True)


if __name__ == '__main__':
    main()

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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosOcpOptions
from time_optimal_simple_bicycle_model import export_time_optimal_simple_bicycle
import numpy as np
from matplotlib import pyplot as plt
import casadi as ca

def create_ddp_opts(Tf, N, M):
    solver_options = AcadosOcpOptions()
    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    solver_options.qp_solver_cond_N = N
    solver_options.hessian_approx = 'GAUSS_NEWTON'
    solver_options.integrator_type = 'ERK'
    solver_options.sim_method_num_steps = M
    solver_options.print_level = 1

    # DDP options
    solver_options.nlp_solver_type = 'DDP'
    solver_options.globalization = 'MERIT_BACKTRACKING'
    solver_options.nlp_solver_max_iter = 200
    solver_options.with_adaptive_levenberg_marquardt = True

    # set prediction horizon
    solver_options.tf = Tf

    return solver_options

def create_sqp_opts(Tf, N, M, globalized:bool = True):
    solver_options = AcadosOcpOptions()
    solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    solver_options.qp_solver_cond_N = N
    solver_options.hessian_approx = 'EXACT'
    solver_options.integrator_type = 'ERK'
    solver_options.sim_method_num_steps = M
    solver_options.print_level = 1

    # DDP options
    solver_options.nlp_solver_type = 'SQP'
    if globalized:
        solver_options.globalization = 'FUNNEL_L1PEN_LINESEARCH'
    else:
        solver_options.globalization = 'FIXED_STEP'
    solver_options.nlp_solver_max_iter = 200
    solver_options.with_adaptive_levenberg_marquardt = False
    solver_options.regularize_method = 'PROJECT'

    # set prediction horizon
    solver_options.tf = Tf

    return solver_options

def solve_acados_ocp(solve_feasibility_problem: bool, acados_opts: AcadosOcpOptions, plotting: bool=False):

    # The flag denotes, if the problem should be transformed into a feasibility
    # problem, or if the unconstrained OCP should be solved.
    SOLVE_FEASIBILITY_PROBLEM = solve_feasibility_problem

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_time_optimal_simple_bicycle()
    ocp.model = model

    Tf = 1.0 # It is a time optimal problem with time scaling
    nx = model.x.rows()
    nu = model.u.rows()
    N = 20
    M = 2 # Needed for integrator

    # set dimensions
    ocp.dims.N = N

    x0 = 0.0
    y0 = 0.0
    theta0 = np.pi/2

    xf = 0.0
    yf = 11.0
    thetaf = np.pi/2

    v_max = 1.0
    delta_max = np.pi/6

    # the 'EXTERNAL' cost type can be used to define general cost terms
    # NOTE: This leads to additional (exact) hessian contributions when using GAUSS_NEWTON hessian.
    if not SOLVE_FEASIBILITY_PROBLEM:
        ocp.cost.cost_type_e = 'EXTERNAL'
        ocp.model.cost_expr_ext_cost_e = model.x[0]
    ###########################################################################
    # Define constraints
    ###########################################################################

    # Initial conditions
    ocp.constraints.lbx_0 = np.array([x0, y0, theta0])
    ocp.constraints.ubx_0 = np.array([x0, y0, theta0])
    ocp.constraints.idxbx_0 = np.array([1,2,3])

    # Path constraints on control
    ocp.constraints.lbu = np.array([-delta_max, 0])
    ocp.constraints.ubu = np.array([+delta_max, +v_max])
    ocp.constraints.idxbu = np.array([0, 1])

    # Terminal constraints
    ocp.constraints.lbx_e = np.array([xf, yf, thetaf])
    ocp.constraints.ubx_e = np.array([xf, yf, thetaf])
    ocp.constraints.idxbx_e = np.array([1,2,3])

    # Convex over Nonlinear Constraints
    # r = ca.SX.sym('r', 1, 1)
    # ocp.model.con_phi_expr = r**2
    # ocp.model.con_r_in_phi = r
    # ocp.model.con_r_expr = (5*(model.x[1]+3)) / (ca.sqrt(1+15*(model.x[2]-5)**2))

    # ocp.constraints.lphi = np.array([0])
    # ocp.constraints.uphi = np.array([1])

    # Otherwise exact Hessian does not exist so far.
    ocp.model.con_h_expr = ((5*(model.x[1]+3)) / (ca.sqrt(1+15*(model.x[2]-5)**2)))**2
    ocp.constraints.lh = np.array([-1])
    ocp.constraints.uh = np.array([1])

    ###########################################################################
    # set solver options
    ###########################################################################
    ocp.solver_options = acados_opts

    if SOLVE_FEASIBILITY_PROBLEM:
        ocp.translate_to_feasibility_problem()

    ocp_solver = AcadosOcpSolver(ocp, json_file = 'hour_glass_acados.json')

    for i in range(N):
        ocp_solver.cost_set(i, "scaling", 1.0)
    sol_X = np.zeros((N+1, nx))
    sol_U = np.zeros((N, nu))

    # Load and set the initial guess
    with open('hour_glass_initial_guess.npy', 'rb') as f:
        X_init = np.load(f)
        U_init = np.load(f)

    for i in range(N):
        ocp_solver.set(i, "x", X_init[:,i])
        ocp_solver.set(i, "u", U_init[:,i])
    ocp_solver.set(N, "x", X_init[:,N])

    # Solve the problem
    status = ocp_solver.solve()

    iter = ocp_solver.get_stats('nlp_iter')

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        sol_X[i,:] = ocp_solver.get(i, "x")
        sol_U[i,:] = ocp_solver.get(i, "u")
    sol_X[N,:] = ocp_solver.get(N, "x")

    if plotting:
        plot_trajectory([X_init, sol_X.T], ["Initial guess", "Solution"])

    return iter


def plot_trajectory(list_X_sol: list, labels: list):

    # Initial states
    x0 = 0
    y0 = 0
    # Final states
    xf = 0
    yf = 11
    # Number of control intervals
    N = 20

    # Plot the x and y solutions + tunnel shape
    plt.figure()
    plt.title('TUNNEL N=' +str(N)+' Start: (' + str("{:.2f}".format(x0)) + ',' + str("{:.2f}".format(y0)) + ') End: (' + str("{:.2f}".format(xf)) + ',' + str("{:.2f}".format(yf)) + ')')
    plt.axis([-8, 5, -1 , 11])

    # Plot the tunnel
    y_co = np.linspace(y0, yf, 1000) # Points along the path
    r = np.sqrt(1+15*(y_co-5)**2)/5    # Distance between the path and the tunnel walls
    plt.plot(-3+r,y_co, 'r-')
    plt.plot(-3-r,y_co, 'r-')
    for i in range(len(y_co)):
        plt.hlines(y= y_co[i], xmin=-3-r[i], xmax=-3+r[i], color='r')
    angle = np.linspace(0,np.pi/4,500, endpoint=True)
    rf1 = np.sqrt(1+15*(yf-5)**2)/5
    rf2 = np.sqrt(1+15*(y0-5)**2)/5
    for i in range(500):
        plt.hlines(y= yf+rf1*np.sin(angle[i]), xmin=-3-rf1*np.cos(angle[i]), xmax=-3+rf1*np.cos(angle[i]), color='r')
        plt.hlines(y= y0-rf2*np.sin(angle[i]), xmin=-3-rf2*np.cos(angle[i]), xmax=-3+rf2*np.cos(angle[i]), color='r')

    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    plt.plot(-3+rf1*cos_angle,yf+rf1*sin_angle, 'r-')
    plt.plot(-3-rf1*cos_angle,yf+rf1*sin_angle, 'r-')
    plt.plot(-3+rf2*cos_angle,y0-rf2*sin_angle, 'r-')
    plt.plot(-3-rf2*cos_angle,y0-rf2*sin_angle, 'r-')

    # Plot x and y coordinates
    for i in range(len(list_X_sol)):
        plt.plot(list_X_sol[i][1,:], list_X_sol[i][2,:], label = labels[i], marker='o', markersize=4)
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.legend()
    plt.show()

def main():
    N = 20
    M = 2 
    Tf = 1.0

    # # Test DDP here
    ddp_opts = create_ddp_opts(Tf, N, M)
    ddp_iter = solve_acados_ocp(True, ddp_opts)
    assert ddp_iter <= 14, "DDP Solver should converge within 14 iterations!"

    # Compare SQP full step with SQP funnel
    sqp_funnel_opts = create_sqp_opts(Tf, N, M, globalized=True)
    sqp_full_step_opts = create_sqp_opts(Tf, N, M, globalized=False)

    sqp_funnel_iter = solve_acados_ocp(False, sqp_funnel_opts)
    sqp_full_step_iter = solve_acados_ocp(False, sqp_full_step_opts)

    assert sqp_full_step_iter == sqp_funnel_iter, "Funnel should only take full steps!"

if __name__ == '__main__':
    main()
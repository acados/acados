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

import sys, json, os
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver, acados_dae_model_json_dump, get_acados_path, get_simulink_default_opts
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum

import matplotlib.pyplot as plt

import casadi as ca

TOL = 1e-7

def main(discretization='shooting_nodes'):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    integrator_type = 'LIFTED_IRK' # ERK, IRK, GNSF, LIFTED_IRK

    if integrator_type == 'GNSF':
        acados_dae_model_json_dump(model)
        # structure detection in MATLAB/Octave -> produces 'pendulum_ode_gnsf_functions.json'
        status = os.system('octave detect_gnsf_from_json.m')
        # load gnsf from json
        with open(model.name + '_gnsf_functions.json', 'r') as f:
            gnsf_dict = json.load(f)
        ocp.gnsf_model = gnsf_dict

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx
    N = 15

    # discretization
    ocp.solver_options.N_horizon = N
    # shooting_nodes = np.linspace(0, Tf, N+1)

    time_steps = np.linspace(0, 1, N+1)[1:]
    time_steps = Tf * time_steps / sum(time_steps)

    shooting_nodes = np.zeros((N+1,))
    for i in range(len(time_steps)):
        shooting_nodes[i+1] = shooting_nodes[i] + time_steps[i]

    # nonuniform discretizations can be defined either by shooting_nodes or time_steps:
    if discretization == 'shooting_nodes':
        ocp.solver_options.shooting_nodes = shooting_nodes
    elif discretization == 'time_steps':
        ocp.solver_options.time_steps = time_steps
    else:
        raise NotImplementedError(f"discretization type {discretization} not supported.")

    # set num_steps
    ocp.solver_options.sim_method_num_steps = 2*np.ones((N,))
    ocp.solver_options.sim_method_num_steps[0] = 3

    # set num_stages
    ocp.solver_options.sim_method_num_stages = 2*np.ones((N,))
    ocp.solver_options.sim_method_num_stages[0] = 4

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])

    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0
    ocp.constraints.idxbu = np.array([0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    # partial condensing settings
    qp_solver_cond_N = 8
    ocp.solver_options.qp_solver_cond_N = qp_solver_cond_N
    ocp.solver_options.qp_solver_cond_block_size = (qp_solver_cond_N) * [1] + [N-((qp_solver_cond_N))]
    # HPIPM settings
    ocp.solver_options.qp_solver_mu0 = 1e3

    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = integrator_type
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.hpipm_mode = "ROBUST"

    ocp.solver_options.globalization = 'FIXED_STEP'

    # set prediction horizon
    ocp.solver_options.tf = Tf

    # Set additional options for Simulink interface:
    simulink_opts = get_simulink_default_opts()
    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json', simulink_opts = simulink_opts, verbose=False)

    simX = np.zeros((N+1, nx))
    simU = np.zeros((N, nu))

    # change options after creating ocp_solver
    ocp_solver.options_set("globalization_fixed_step_length", 0.99999)
    ocp_solver.options_set("tol_eq", TOL)
    ocp_solver.options_set("tol_stat", TOL)
    ocp_solver.options_set("tol_ineq", TOL)
    ocp_solver.options_set("tol_comp", TOL)

    # initialize solver
    for i in range(N):
        ocp_solver.set(i, "x", x0)
    status = ocp_solver.solve()

    # get primal solution
    for i in range(N):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N,:] = ocp_solver.get(N, "x")

    # get condensed Hessian
    pcond_H = []
    pcond_Q = []
    pcond_R = []
    pcond_S = []
    for i in range(ocp.solver_options.qp_solver_cond_N+1):
        pcond_Q.append(ocp_solver.get_from_qp_in(i, "pcond_Q"))
        pcond_R.append(ocp_solver.get_from_qp_in(i, "pcond_R"))
        pcond_S.append(ocp_solver.get_from_qp_in(i, "pcond_S"))

        pcond_RSQ = ca.blockcat(pcond_Q[-1], pcond_S[-1].T, pcond_S[-1], pcond_R[-1]).full()
        # copy lower triangular part to upper triangular part
        pcond_RSQ = np.tril(pcond_RSQ) + np.tril(pcond_RSQ, -1).T
        pcond_H.append(pcond_RSQ)
        eig_values, _ = np.linalg.eig(pcond_RSQ)
        print(f"eigen values of partially condenesed Hessian block {i}: {eig_values}")

    # pcond_H = ocp_solver.get_from_qp_in(ocp.solver_options.qp_solver_cond_N, "pcond_H")
    # print("pcond_H", pcond_H)
    # pcond_H_mat = scipy.linalg.block_diag(*pcond_H)
    # plt.spy(pcond_H_mat)
    # plt.show()

    # check dimensions of partially condensed matrices
    assert pcond_Q[0].shape == (0, 0)
    for i in range(1, ocp.solver_options.qp_solver_cond_N+1):
        assert pcond_Q[i].shape == (nx, nx)
    for i in range(ocp.solver_options.qp_solver_cond_N+1):
        block_size = ocp.solver_options.qp_solver_cond_block_size[i]
        assert pcond_R[i].shape == (nu*block_size, nu*block_size)

    print("inequality multipliers at stage 1")
    print(ocp_solver.get(1, "lam")) # inequality multipliers at stage 1
    print("multipliers of dynamic conditions between stage 1 and 2")
    print(ocp_solver.get(1, "pi")) # multipliers of dynamic conditions between stage 1 and 2

    # initialize ineq multipliers at stage 1
    ocp_solver.set(1, "lam", np.zeros(2,))

    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")

    # timings
    time_tot = ocp_solver.get_stats("time_tot")
    time_lin = ocp_solver.get_stats("time_lin")
    time_sim = ocp_solver.get_stats("time_sim")
    time_qp = ocp_solver.get_stats("time_qp")

    print(f"timings OCP solver: total: {1e3*time_tot}ms, lin: {1e3*time_lin}ms, sim: {1e3*time_sim}ms, qp: {1e3*time_qp}ms")
    # print("simU", simU)
    # print("simX", simX)
    iterate_filename = f'final_iterate_{discretization}.json'
    ocp_solver.store_iterate(filename=iterate_filename, overwrite=True)

    if status not in [0, 2]:
        ocp_solver.store_iterate()
        ocp_solver.dump_last_qp_to_json()
        ocp_solver.print_statistics()
        raise Exception(f'acados returned status {status}.')

    plot_pendulum(shooting_nodes, Fmax, simU, simX, latexify=False)
    del ocp_solver


if __name__ == "__main__":
    discretizations = ['shooting_nodes', 'time_steps']
    for discretization in discretizations:
        main(discretization=discretization)

    import json
    # compare iterates
    iterate_filename = f'final_iterate_shooting_nodes.json'
    with open(iterate_filename, 'r') as f:
        iterate0 = json.load(f)

    iterate_filename = f'final_iterate_time_steps.json'
    with open(iterate_filename, 'r') as f:
        iterate1 = json.load(f)

    assert iterate1.keys() == iterate0.keys()

    for k in iterate0:
        assert(len(iterate0[k]) == len(iterate1[k]))
        if len(iterate0[k]) > 0:
            diff = np.max(np.abs(np.array(iterate0[k]) - np.array(iterate1[k])))
            if diff > TOL:
                raise Exception(f'results for {k}, do not match up to accuracy {TOL}, diff is {diff}')

    print(f'Check passed: discretization formulations are equivalent: {discretizations}')
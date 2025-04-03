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
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver, acados_dae_model_json_dump, get_acados_path, get_simulink_default_opts
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum

import matplotlib.pyplot as plt

import casadi as ca

TOL = 1e-7

def main(nlp_solver_type="SQP"):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    integrator_type = 'ERK'

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx
    N = 15

    # discretization
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

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

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'

    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = integrator_type
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = nlp_solver_type

    ocp.solver_options.rti_log_residuals = 1
    ocp.solver_options.rti_log_only_available_residuals = 1

    # create ocp solver
    solver = AcadosOcpSolver(ocp)

    if nlp_solver_type == "SQP":
        solver.solve()
        solver.print_statistics()
        res_stat_all = solver.get_stats("res_stat_all")
        res_eq_all = solver.get_stats("res_eq_all")
        res_ineq_all = solver.get_stats("res_ineq_all")
        res_comp_all = solver.get_stats("res_comp_all")
        nlp_iter = solver.get_stats('nlp_iter')

        res_init = solver.get_initial_residuals()
        print(f"residuals: {res_init}")
        assert res_init[0] == res_stat_all[0]
        assert res_init[1] == res_eq_all[0]
        assert res_init[2] == res_ineq_all[0]
        assert res_init[3] == res_comp_all[0]

    elif nlp_solver_type == "SQP_RTI":
        TOL = 1e-6
        res_stat_all = []
        res_eq_all = []
        res_ineq_all = []
        res_comp_all = []
        for nlp_iter in range(20):
            solver.solve()
            solver.print_statistics()
            res_stat_all.append(solver.get_stats("res_stat_all")[0])
            res_eq_all.append(solver.get_stats("res_eq_all")[0])
            res_ineq_all.append(solver.get_stats("res_ineq_all")[0])
            res_comp_all.append(solver.get_stats("res_comp_all")[0])
            #
            # quick way to get all residuals
            res_all = solver.get_initial_residuals()
            print(f"res_all: {res_all}")

            assert res_all[0] == res_stat_all[-1]
            assert res_all[1] == res_eq_all[-1]
            assert res_all[2] == res_ineq_all[-1]
            assert res_all[3] == res_comp_all[-1]
            if nlp_iter > 0:
                # veryfiy that the residuals were "predicted" correctly
                assert res_next[0] == res_stat_all[-1]
                assert res_next[1] == res_eq_all[-1]
                assert res_next[2] == res_ineq_all[-1]
                assert res_next[3] == res_comp_all[-1]
            # overwrite res_next
            res_next = solver.get_residuals()
            print(f"res_next: {res_next}")

            if max([res_stat_all[-1], res_eq_all[-1], res_ineq_all[-1], res_comp_all[-1]]) < TOL:
                break

    print(f"obtained residuals after {nlp_iter} iterations")
    for i in range(nlp_iter):
        print(f"{i}: {res_stat_all[i]:e} {res_eq_all[i]:e} {res_ineq_all[i]:e} {res_comp_all[i]:e}")

    del solver
    return res_stat_all, res_eq_all, res_ineq_all, res_comp_all



if __name__ == "__main__":
    res_stat_all_rti, res_eq_all_rti, res_ineq_all_rti, res_comp_all_rti = main(nlp_solver_type="SQP_RTI")
    res_stat_all_sqp, res_eq_all_sqp, res_ineq_all_sqp, res_comp_all_sqp = main(nlp_solver_type="SQP")

    # check that number of iterations coincide
    n_iter_rti = len(res_stat_all_rti)
    n_iter_sqp = len(res_stat_all_sqp)
    if n_iter_rti != n_iter_sqp:
        raise Exception(f"number of iterations differ: rti {n_iter_rti}, sqp {n_iter_sqp}")

    # check that residuals conincide
    for name, r_rti, r_sqp in [("stat", res_stat_all_rti, res_stat_all_sqp),
                         ("eq", res_eq_all_rti, res_eq_all_sqp),
                         ("ineq", res_ineq_all_rti, res_ineq_all_sqp),
                         ("comp", res_comp_all_rti, res_comp_all_sqp)]:
        if not np.allclose(r_rti, r_sqp, atol=1e-6):
            print(f"rti: {r_rti}")
            print(f"sqp: {r_sqp}")
            raise Exception(f"residuals {name} do not coincide")

    print("SUCCESS: residuals of SQP and SQP_RTI coincide")

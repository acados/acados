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
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import vertcat
from test_sens_forw_p import export_pendulum_model_with_M_param

def main():
    model = export_pendulum_model_with_M_param()
    model_name_base = model.name

    # OCP Dimensions
    N = 10
    Tf = 1.0
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    np_param = model.p.size()[0]

    for integrator_type in ['ERK', 'IRK']:
        print(f"\n--- Testing OCP with Integrator: {integrator_type} ---")

        ocp = AcadosOcp()
        ocp.model = model
        ocp.model.name = f"{model_name_base}_ocp_{integrator_type}"
        ocp.solver_options.N_horizon = N
        ocp.solver_options.tf = Tf
        ocp.parameter_values = np.array([1.0])

        # Cost
        Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
        R_mat = 2 * np.diag([1e-2])

        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.model.cost_y_expr = vertcat(model.x, model.u)
        ocp.cost.yref = np.zeros((nx+nu,))
        ocp.cost.W = np.diag(np.concatenate([np.diag(Q_mat), np.diag(R_mat)]))

        ocp.cost.cost_type_e = 'NONLINEAR_LS'
        ocp.model.cost_y_expr_e = model.x
        ocp.cost.yref_e = np.zeros((nx,))
        ocp.cost.W_e = Q_mat

        # Constraints
        Fmax = 80
        ocp.constraints.lbu = np.array([-Fmax])
        ocp.constraints.ubu = np.array([+Fmax])
        ocp.constraints.idxbu = np.array([0])
        ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

        # Solver Options
        ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
        ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
        ocp.solver_options.integrator_type = integrator_type
        ocp.solver_options.nlp_solver_type = 'SQP'
        ocp.solver_options.sens_forw_p = True
        ocp.solver_options.sim_method_num_stages = 4
        ocp.solver_options.sim_method_num_steps = 10

        ocp_solver = AcadosOcpSolver(ocp)

        # Set Parameters
        p_val = np.array([1.0])
        for i in range(N + 1):
            ocp_solver.set(i, "p", p_val)

        if ocp_solver.solve() != 0:
            raise Exception('acados returned error status.')

        # -------------------------------------------------------------------------
        # Verification: Compare OCP S_p at stage 0 against Integrator FD
        # -------------------------------------------------------------------------

        # Verify at stage 1 (result of first shooting interval)
        stage = 0
        S_p_ocp = ocp_solver.get(stage, "S_p")

        # Get operating point at stage 0
        x0_op = ocp_solver.get(0, "x")
        u0_op = ocp_solver.get(0, "u")

        # Setup Sim for ground truth
        sim = AcadosSim()
        sim.model = model
        # Use a unique name for sim as well to be safe
        sim.model.name = f"{model_name_base}_sim_{integrator_type}"
        sim.solver_options.T = Tf / N
        sim.solver_options.integrator_type = integrator_type
        sim.solver_options.num_stages = 4
        sim.solver_options.num_steps = 10

        sim.parameter_values = p_val
        sim_solver = AcadosSimSolver(sim)

        # Finite Differences on Sim
        FD_epsilon = 1e-6
        S_p_fd = np.zeros((nx, np_param))

        sim_solver.set('x', x0_op)
        sim_solver.set('u', u0_op)
        sim_solver.set('p', p_val)
        sim_solver.solve()
        x_next_nom = sim_solver.get('x')

        for j in range(np_param):
            p_pert = p_val.copy()
            p_pert[j] += FD_epsilon

            sim_solver.set('x', x0_op)
            sim_solver.set('u', u0_op)
            sim_solver.set('p', p_pert)
            sim_solver.solve()
            x_next_pert = sim_solver.get('x')

            S_p_fd[:, j] = (x_next_pert - x_next_nom) / FD_epsilon

        error = np.max(np.abs(S_p_ocp - S_p_fd))
        print(f"Stage {stage} Max Error ({integrator_type} OCP S_p vs Sim FD): {error:.2e}")

        # print("OCP S_p:",  S_p_ocp)
        # print("Sim FD:",  S_p_fd)

        if error > 1e-5:
            raise Exception(f"Failure: OCP parameter sensitivities do not match dynamics for {integrator_type}.")

        del ocp_solver
        del sim_solver

    print("\nSuccess: OCP parameter sensitivities verified for both ERK and IRK.")

if __name__ == '__main__':
    main()
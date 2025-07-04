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

# reference : "Towards Time-optimal Tunnel-following for Quadrotors", Jon Arrizabalaga et al.

import casadi as ca
from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, AcadosSimSolver, ACADOS_INFTY
import scipy.linalg

from common import *
from sys_dynamics import SysDyn

class AcadosCustomOcp:

    def __init__(self):
        self.nx = 0
        self.nu = 0
        self.ny = 0
        self.ns = 0

        self.ocp = None,
        self.solver = None,
        self.integrator = None
        self.sysModel = None

        self.zeta_0 = None
        self.zeta_N = None
        self.u_N = None


    def setup_acados_ocp(self):
        '''Formulate acados OCP'''

        # create casadi symbolic expressions
        sysModel = SysDyn()
        self.sysModel = sysModel

        zeta_f, dyn_f, u, proj_constr, dyn_fn = sysModel.SetupOde()
        self.zeta_0 = np.copy(init_zeta)

        # create Acados model
        ocp = AcadosOcp()
        model_ac = AcadosModel()
        model_ac.f_expl_expr = dyn_f
        model_ac.x = zeta_f
        model_ac.u = u
        model_ac.name = "drone_FrenSer"
        ocp.model = model_ac

        # set dimensions
        ocp.solver_options.N_horizon = N
        self.nx = model_ac.x.size()[0]
        self.nu = model_ac.u.size()[0]
        ny = self.nx + self.nu

        self.zeta_N = ca.repmat(np.reshape(self.zeta_0, (self.nx,1)), 1, N+1)
        self.u_N = ca.repmat(U_REF, 1, N)

        # continuity constraints
        ocp.constraints.x0  = self.zeta_0

        # formulate cost function
        ocp.cost.cost_type = "NONLINEAR_LS"
        ocp.model.cost_y_expr = ca.vertcat(model_ac.x, model_ac.u)
        ocp.cost.yref = np.array([ 0.2, 0, 0,
                                   1, 0, 0, 0,
                                   0, 0, 0,
                                   0, 0, 0,
                                   0, 0, 0,
                                  U_HOV, U_HOV, U_HOV, U_HOV,
                                  0, 0, 0, 0])
        ocp.cost.W = scipy.linalg.block_diag(Q, R)

        ocp.cost.cost_type_e = "NONLINEAR_LS"
        ocp.model.cost_y_expr_e = model_ac.x
        ocp.cost.yref_e = np.array([ 0.2, 0, 0,
                                     1, 0, 0, 0,
                                     0, 0, 0,
                                     0, 0, 0,
                                     0, 0, 0,
                                     U_HOV, U_HOV, U_HOV, U_HOV])
        ocp.cost.W_e = Qn

        # formulate inquality constraints

        # constrain AGV dynamics : acceleration, angular velocity (convex ?, Non-linear)
        dyn_constr_eqn = []
        dyn_constr_eqn = ca.vertcat(dyn_constr_eqn , proj_constr)

        ineq_constr_eqn = []
        ineq_constr_eqn = ca.vertcat(ineq_constr_eqn, dyn_constr_eqn)

        model_ac.con_h_expr = ineq_constr_eqn
        model_ac.con_h_expr_e = ineq_constr_eqn

        # inequality bounds
        nh = model_ac.con_h_expr.shape[0]

        # constrain controls
        # lbu = [0] * self.nu;      ubu = [0] * self.nu

        # # Control bounds ( Affects horizon quality before switch)
        # lbu[0] = OHM_MIN;         ubu[0] = OHM_MAX
        # lbu[1] = OHM_MIN;         ubu[1] = OHM_MAX

        # ocp.constraints.lbu = np.array(lbu)
        # ocp.constraints.ubu = np.array(ubu)
        # ocp.constraints.idxbu = np.array([0, 1])

        # Bounds on path constraints (inequality)
        lh = np.zeros(nh);        uh = np.zeros(nh)
        lh[:] = -ACADOS_INFTY;             uh[:] = 1

        ocp.constraints.lh = lh
        ocp.constraints.uh = uh

        ocp.constraints.lh_e = lh
        ocp.constraints.uh_e = uh

        # configure itegrator and QP solver
        ocp.solver_options.integrator_type = "ERK"
        ocp.solver_options.tf = Tf
        ocp.solver_options.sim_method_num_stages = 4
        ocp.solver_options.sim_method_num_steps = 1
        # ocp.solver_options.collocation_type = 'GAUSS_RADAU_IIA'
        # ocp.solver_options.time_steps = time_steps
        # ocp.solver_options.shooting_nodes = shooting_nodes

        ocp.solver_options.qp_solver = "PARTIAL_CONDENSING_HPIPM"#"PARTIAL_CONDENSING_HPIPM" #"FULL_CONDENSING_HPIPM" #"PARTIAL_CONDENSING_HPIPM"
        ocp.solver_options.hessian_approx =  "GAUSS_NEWTON"#"EXACT",
        # ocp.solver_options.cost_discretization ="INTEGRATOR"
        ocp.solver_options.qp_solver_cond_N = int(N/2)
        ocp.solver_options.nlp_solver_type = "SQP_RTI"
        ocp.solver_options.tol = 1e-3

        # create solver
        self.ocp = ocp
        self.solver = AcadosOcpSolver(ocp, json_file = "planner_ocp.json")
        self.integrator = AcadosSimSolver(ocp)

        return True


    def solve_and_sim(self):
        '''Solve the OCP with multiple shooting, and forward simulate with RK4'''

        # Integrate ODE model to get CL estimate (no measurement noise)
        u_0 = self.solver.solve_for_x0(self.zeta_0)
        self.zeta_0 = self.integrator.simulate(x=self.zeta_0, u=u_0)

        u_0 = self.solver.solve_for_x0(self.zeta_0)
        self.zeta_N = np.reshape(self.solver.get(0, "x"), (self.nx, 1))
        for i in range(1, N +1):
            zeta_i = np.reshape(self.solver.get(i, "x"), (self.nx, 1))
            self.zeta_N = np.concatenate((self.zeta_N, zeta_i), axis = 1)

        self.u_N[:, 0] = u_0


    def cost_update_ref(self, zeta_0, u_ref):

        s0 = float(zeta_0[0])
        if s0 >= S_MAX:
            return True

        sref =  s0 + S_REF
        srefDot = S_REF / Tf
        for j in range(N):
            sref_j = s0 + (sref - s0) * j /N
            yref = np.array([sref_j, 0, 0, 0, 0, 0, 0, srefDot, 0, 0, 0, 0, 0, 0, 0, 0, U_HOV, U_HOV, U_HOV, U_HOV, 0, 0, 0, 0])
            self.solver.set(j, "yref", yref)

        yref = np.array([sref_j, 0, 0, 0, 0, 0, 0, srefDot, 0, 0, 0, 0, 0, 0, 0, 0, U_HOV, U_HOV, U_HOV, U_HOV])
        self.solver.set(N, "yref", yref)


    def get_cost(self):
        cost = self.solver.get_cost()
        return cost


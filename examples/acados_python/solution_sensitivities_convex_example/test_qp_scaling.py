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
from acados_template import AcadosOcp, AcadosModel, AcadosOcpSolver, latexify_plot
latexify_plot()

def export_parametric_ocp() -> AcadosOcp:

    model = AcadosModel()
    model.x = ca.SX.sym("x", 1)
    model.p_global = ca.SX.sym("p_global", 1)
    model.disc_dyn_expr = model.x
    model.cost_expr_ext_cost = (model.x - model.p_global**2)**2
    model.cost_expr_ext_cost_e = 0
    model.name = "non_ocp"
    ocp = AcadosOcp()
    ocp.model = model

    ocp.constraints.lbx_0 = np.array([-1.0])
    ocp.constraints.ubx_0 = np.array([1.0])
    ocp.constraints.idxbx_0 = np.array([0])

    ocp.cost.cost_type = "EXTERNAL"
    ocp.solver_options.integrator_type = "DISCRETE"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.N_horizon = 1
    ocp.solver_options.tf = 1.0

    ocp.p_global_values = np.zeros((1,))

    return ocp

def main():
    p_test = [-1.5]

    ocp = export_parametric_ocp()
    ocp.solver_options.qp_scaling_type = "OBJECTIVE_GERSHGORIN"
    ocp.solver_options.nlp_solver_max_iter = 2

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)


    for i, p in enumerate(p_test):
        p_val = np.array([p])

        ocp_solver.set_p_global_and_precompute_dependencies(p_val)
        status = ocp_solver.solve()
        iterate = ocp_solver.store_iterate_to_flat_obj()
        print(f"iterate: {iterate.x}, {iterate.u}, {iterate.lam}, {iterate.pi}")

        if status != 0:
            # ocp_solver.print_statistics()
            raise Exception(f"OCP solver returned status {status} at {i}th p value {p}.")



if __name__ == "__main__":
    main()

#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import os
import numpy as np
import casadi as ca
from acados_template import AcadosOcp, AcadosModel, AcadosOcpSolver, ocp_get_default_cmake_builder
from subprocess import call

def export_parametric_nlp() -> AcadosOcp:

    model = AcadosModel()
    model.x = ca.SX.sym("x", 1)
    model.p_global = ca.SX.sym("p_global", 1)
    model.cost_expr_ext_cost_e = (model.x - model.p_global**2)**2
    model.name = "non_ocp"
    ocp = AcadosOcp()
    ocp.model = model

    ocp.constraints.lbx_e = np.array([-1.0])
    ocp.constraints.ubx_e = np.array([1.0])
    ocp.constraints.idxbx_e = np.array([0])

    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.solver_options.qp_solver = "FULL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.N_horizon = 0

    ocp.p_global_values = np.zeros((1,))
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    return ocp


def main():
    ocp = export_parametric_nlp()

    ocp.code_gen_options.json_file = 'test_nlp.json'
    cmake_builder = ocp_get_default_cmake_builder()
    cmake_builder.options_on = ['BUILD_EXAMPLE_STATIC', 'BUILD_ACADOS_OCP_SOLVER_STATIC_LIB']
    AcadosOcpSolver.generate(ocp, cmake_builder=cmake_builder)
    code_export_dir = ocp.code_gen_options.code_export_directory
    cmake_builder.exec(code_export_dir, verbose=True)

    call(os.path.join(os.getcwd(), ocp.code_export_directory, f'main_{ocp.name}'))


if __name__ == "__main__":
    main()

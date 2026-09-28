#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

import numpy as np
from acados_template import AcadosOcpSolver
from non_ocp_example import create_parametric_nlp

# test barrier QP residuals
def test_barrier_qp_residual():
    p_nominal = 0.0
    np_test = 50
    p_test = np.linspace(p_nominal, p_nominal + 2, np_test)

    ocp = create_parametric_nlp()
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration
    # test doesnt need solution sensitivities
    ocp.code_gen_options.with_solution_sens_wrt_params = False
    ocp.code_gen_options.with_value_sens_wrt_params = False

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    for tau in [0.0, 1e-2, 1e-3]:
        ocp_solver.options_set("tau_min", tau)
        for i, p in enumerate(p_test):
            p_val = np.array([p])

            ocp_solver.set_p_global_and_precompute_dependencies(p_val)
            status = ocp_solver.solve()
            # ocp_solver.print_statistics()
            if status != 0:
                raise Exception(f"OCP solver returned status {status} at {i}th p value {p}, {tau=}.")
            qp_residuals = ocp_solver.get_stats("qp_residuals")
            if any(qp_residuals > ocp.solver_options.tol):
                raise Exception(f"QP residuals too high: {qp_residuals} at {i}th p value {p}, {tau=}.")
        print(f"test_barrier_qp_residual passed for tau={tau}")
        print('last solver stats:')
        ocp_solver.print_statistics()



if __name__ == "__main__":
    test_barrier_qp_residual()

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
from acados_template import AcadosOcpSolver
from non_ocp_example import export_parametric_nlp

# test barrier QP residuals
def test_barrier_qp_residual():
    p_nominal = 0.0
    np_test = 50
    p_test = np.linspace(p_nominal, p_nominal + 2, np_test)

    ocp = export_parametric_nlp()
    ocp.solver_options.qp_solver_t0_init = 0
    ocp.solver_options.nlp_solver_ext_qp_res = 1
    ocp.solver_options.nlp_solver_max_iter = 2 # QP should converge in one iteration
    # test doesnt need solution sensitivities
    ocp.solver_options.with_solution_sens_wrt_params = False
    ocp.solver_options.with_value_sens_wrt_params = False

    ocp_solver = AcadosOcpSolver(ocp, json_file="parameter_augmented_acados_ocp.json", verbose=False)

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

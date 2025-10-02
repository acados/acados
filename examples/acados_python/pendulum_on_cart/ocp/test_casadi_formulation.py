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

import sys
sys.path.insert(0, '../common')

import numpy as np
import casadi as ca
from typing import Union

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from ocp_example_cost_formulations import formulate_ocp, T_HORIZON, FMAX

from utils import plot_pendulum

def raise_test_failure_message(msg: str):
    # print(f"ERROR: {msg}")
    raise Exception(msg)

def main(cost_version="LS", constraint_version='h', casadi_solver_name="ipopt", use_acados_hessian=False):
    ocp = formulate_ocp(cost_version=cost_version, constraint_version=constraint_version)
    ocp.solver_options.tf = T_HORIZON

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    # initialize solver
    ocp_solver.load_iterate_from_obj(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    ocp_solver.print_statistics()
    if status != 0:
        raise_test_failure_message(f"acados solver returned status {status}.")
    # get solution
    result_acados = ocp_solver.store_iterate_to_obj()
    result_acados_flat = ocp_solver.store_iterate_to_flat_obj()

    ## solve using casadi
    casadi_solver_opts = {}
    if casadi_solver_name == "fatrop":
        casadi_solver_opts["expand"] = True
        casadi_solver_opts["fatrop"] = {"mu_init": 0.1}
        casadi_solver_opts["structure_detection"] = "auto"
        casadi_solver_opts["debug"] = True

    if use_acados_hessian:
        # ocp.solver_options.hessian_approx = "EXACT"
        ocp.solver_options.fixed_hess = False

    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False, solver=casadi_solver_name, casadi_solver_opts=casadi_solver_opts, use_acados_hessian=use_acados_hessian)
    casadi_ocp_solver.load_iterate_from_obj(result_acados)
    status = casadi_ocp_solver.solve()
    print(f"casadi solver returned status {status}.")
    result_casadi = casadi_ocp_solver.store_iterate_to_obj()
    result_casadi_flat = casadi_ocp_solver.store_iterate_to_flat_obj()

    nlp_iter_ca = casadi_ocp_solver.get_stats("nlp_iter")
    print(f"Casadi solver finished after {nlp_iter_ca} iterations.")

    # evaluate difference
    # TODO: set solver tolerance and reduce it here.
    test_tol = 5e-5
    for field in ["x", "u", "pi", "lam"]:
        diff = np.linalg.norm(getattr(result_acados_flat, field) - getattr(result_casadi_flat, field))
        print(f"Difference between acados and casadi solution for {field:4}: {diff:.3e}")
        if diff > test_tol:
            raise_test_failure_message(f"Test failed: difference between acados and casadi solution for {field} should be smaller than {test_tol}, but got {diff}.")

    # additional checks:
    if cost_version == "LS" and constraint_version == "h":
        if use_acados_hessian:
            if nlp_iter_ca < 10:
                raise_test_failure_message(f"Expected more iterations for casadi solver with acados (GN) Hessian, but got {nlp_iter_ca} iterations.")
        else:
            if nlp_iter_ca > 10:
                raise_test_failure_message(f"Expected less iterations for casadi solver with Hessian, but got {nlp_iter_ca} iterations.")

    # plot_pendulum(ocp.solver_options.shooting_nodes, FMAX,
    #               np.array(result_casadi.u_traj), np.array(result_casadi.x_traj), latexify=False)
    del ocp_solver, casadi_ocp_solver


if __name__ == "__main__":
    # TODO: refactor to create acados solver only once fore each fomulation version.
    main(cost_version="LS", constraint_version='bu')
    main(cost_version="LS", constraint_version='h', casadi_solver_name="ipopt", use_acados_hessian=True)
    main(cost_version="LS", constraint_version='h', casadi_solver_name="ipopt", use_acados_hessian=False)
    main(cost_version="NLS", constraint_version='h', casadi_solver_name="ipopt", use_acados_hessian=True)
    main(cost_version="NLS", constraint_version='h', casadi_solver_name="ipopt", use_acados_hessian=False)
    main(cost_version="LS", constraint_version='bu', casadi_solver_name="fatrop")
    main(cost_version="LS", constraint_version='h', casadi_solver_name="fatrop")

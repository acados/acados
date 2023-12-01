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


from acados_template import AcadosOcpSolver, AcadosOcp
from piecewiese_polynomial_control_example import create_ocp_solver, create_mocp_solver

import numpy as np

def main():
    # settings:
    cost_type='NONLINEAR_LS'
    explicit_symmetric_penalties=True
    penalty_type='L2'
    N_horizon = 20
    N_1 = 2
    degrees_u_polynom = [0, 1]
    nlp_solver_max_iter = 1

    ocp_solver_1: AcadosOcpSolver
    ocp_solver_1, evaluate_polynomial_u_fun = create_ocp_solver(cost_type, N_horizon, degrees_u_polynom[0], explicit_symmetric_penalties=explicit_symmetric_penalties, penalty_type=penalty_type, nlp_solver_max_iter=nlp_solver_max_iter)

    ocp_solver_2: AcadosOcpSolver
    ocp_solver_2, evaluate_polynomial_u_fun = create_ocp_solver(cost_type, N_horizon, degrees_u_polynom[1], explicit_symmetric_penalties=explicit_symmetric_penalties, penalty_type=penalty_type, nlp_solver_max_iter=nlp_solver_max_iter)

    mocp_solver: AcadosOcpSolver
    mocp_solver, _ = create_mocp_solver(cost_type, [N_1, N_horizon-N_1], degrees_u_polynom, explicit_symmetric_penalties=explicit_symmetric_penalties, penalty_type=penalty_type, nlp_solver_max_iter=nlp_solver_max_iter)

    # call solvers
    # mocp_solver.store_iterate('iter_init_mocp.json', overwrite=True)
    for solver in [ocp_solver_1, ocp_solver_2, mocp_solver]:
        solver.solve()
        solver.print_statistics()
    # mocp_solver.store_iterate('iter_1_mocp.json', overwrite=True)
    # ocp_solver_1.store_iterate('iter_1_ocp1.json', overwrite=True)
    # ocp_solver_2.store_iterate('iter_1_ocp2.json', overwrite=True)

    if len(ocp_solver_1.get(0, 'u')) != degrees_u_polynom[0]+1:
        raise Exception("ocp_solver_1: returned u has wrong dimension")
    if len(ocp_solver_2.get(0, 'u')) != degrees_u_polynom[1]+1:
        raise Exception("ocp_solver_2: returned u has wrong dimension")

    # compare QPs
    print("Comparing QP: Phase 1")
    compare_qp_fields(mocp_solver, ocp_solver_1, range(N_1))
    print("Comparing QP: Phase 2")
    compare_qp_fields(mocp_solver, ocp_solver_2, range(N_1, N_horizon))


def compare_qp_fields(mocp_solver: AcadosOcpSolver, ocp_solver: AcadosOcpSolver, indices):
    for i in indices:
        for field in ['A', 'B', 'b', 'Q', 'R', 'S', 'q', 'r']: #, 'C', 'D']: #mocp_solver.__qp_dynamics_fields:
            val_mocp = mocp_solver.get_from_qp_in(i, field)
            val_ocp = ocp_solver.get_from_qp_in(i, field)
            max_diff = np.max(np.abs(val_mocp - val_ocp))
            if not np.allclose(val_mocp, val_ocp, atol=1e-6):
                print(f"\nfield {field} at stage {i} differs with {max_diff=}, \n {val_mocp=} \n {val_ocp=}\n")
                raise Exception(f"field {field} at stage {i} differs with {max_diff=}, \n {val_mocp=} \n {val_ocp=}\n")

if __name__ == "__main__":
    main()

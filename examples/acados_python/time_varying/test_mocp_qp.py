#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


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

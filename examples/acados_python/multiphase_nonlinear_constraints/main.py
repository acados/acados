#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

from acados_template import AcadosOcpSolver
from create_mocp import create_mocp


def main(soften_h=False, qp_solver='FULL_CONDENSING_QPOASES'):
    multiphase_ocp = create_mocp(soften_h=soften_h, qp_solver=qp_solver)

    ocp_solver = AcadosOcpSolver(multiphase_ocp)

    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    ocp_solver.reset(reset_qp_solver_mem=True, reset_numerical_values=True, reset_solver_options=True, reset_x_to_x0_bar=True)
    assert status == 0, f'acados returned status {status}'

if __name__ == '__main__':
    main()
    main(qp_solver="FULL_CONDENSING_HPIPM", soften_h=True)

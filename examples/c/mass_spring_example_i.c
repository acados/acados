/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// external
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// acados
#include <acados_c/ocp_qp.h>
#include <acados_c/options.h>
#include <acados_c/legacy_create.h>
#include <acados/utils/print.h>

#ifndef ACADOS_WITH_QPDUNES
#define ELIMINATE_X0
#endif
#define GENERAL_CONSTRAINT_AT_TERMINAL_STAGE

#define NREP 100

#include "./mass_spring.c"

int main() {

    printf("\n\n\n mass spring example: acados ocp_qp solvers\n\n\n");

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_in->dim);

    // choose ocp qp solvers
    ocp_qp_solver_t ocp_qp_solvers[] = {
        PARTIAL_CONDENSING_HPIPM,
        FULL_CONDENSING_HPIPM,
        FULL_CONDENSING_QPOASES,
    };

    // choose values for N2 in partial condensing solvers
    int num_N2_values = 3;
    int N2_values[3] = {15, 5, 3};

    for (int ii = 0; ii < 3; ii++) {
        ocp_qp_solver_plan plan;
        plan.qp_solver = ocp_qp_solvers[ii];

        void *args = ocp_qp_create_args(&plan, qp_in->dim);

        for (int jj = 0; jj < num_N2_values; jj++)
        {
            int N2 = N2_values[jj];

            switch (plan.qp_solver) {
                case PARTIAL_CONDENSING_HPIPM:
                    printf("\nPartial condensing + HPIPM (N2 = %d):\n\n", N2);
                    set_option_int(args, "sparse_hpipm.N2", N2);
                    set_option_int(args, "sparse_hpipm.max_iter", 30);
                    break;
                case FULL_CONDENSING_HPIPM:
                    printf("\nFull condensing + HPIPM:\n\n");
                    break;
                case FULL_CONDENSING_QPOASES:
                    printf("\nFull condensing + QPOASES:\n\n");
                    break;
                default:
                    printf("Solver not available\n");
                    exit(1);
                    break;
            }

            ocp_qp_solver *qp_solver = ocp_qp_create(&plan, qp_in->dim, args);

            int acados_return = 0;

            ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
            ocp_qp_info min_info = {ACADOS_POS_INFTY, ACADOS_POS_INFTY, ACADOS_POS_INFTY, ACADOS_POS_INFTY, INT_MAX};

            // run QP solver NREP times and record min timings
            for (int rep = 0; rep < NREP; rep++)
            {
                acados_return += ocp_qp_solve(qp_solver, qp_in, qp_out);

                min_info.num_iter = fmin(min_info.num_iter, info->num_iter);
                min_info.total_time = fmin(min_info.total_time, info->total_time);
                min_info.condensing_time = fmin(min_info.condensing_time, info->condensing_time);
                min_info.solve_QP_time = fmin(min_info.solve_QP_time, info->solve_QP_time);
                min_info.interface_time = fmin(min_info.interface_time, info->interface_time);
            }

            // Compute infinity norm of residuals
            ocp_qp_res *qp_res = create_ocp_qp_res(qp_in->dim);
            ocp_qp_res_ws *res_ws = create_ocp_qp_res_ws(qp_in->dim);
            compute_ocp_qp_res(qp_in, qp_out, qp_res, res_ws);

            double res[4], max_res = 0.0;
            compute_ocp_qp_res_nrm_inf(qp_res, res);
            for (int ii = 0; ii < 4; ii++)
                max_res = fmax(res[ii], max_res);

            printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

            print_ocp_qp_info(&min_info);

            free(qp_solver);
        }
        free(args);
    }

    free(qp_in);
    free(qp_out);

    return 0;
}

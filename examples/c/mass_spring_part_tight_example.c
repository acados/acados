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
#include <stdio.h>
#include <stdlib.h>
// acados
#include <acados_c/ocp_qp.h>
#include <acados_c/options.h>
#include <acados_c/legacy_create.h>
// NOTE(nielsvd): required to cast memory etc. should go.
#include <acados/utils/print.h>
#include <acados/ocp_qp/ocp_qp_sparse_solver.h>
#include <acados/ocp_qp/ocp_qp_full_condensing_solver.h>
#include <acados/ocp_qp/ocp_qp_hpipm.h>
#ifdef ACADOS_WITH_HPMPC
#include <acados/ocp_qp/ocp_qp_hpmpc.h>
#endif
#ifdef ACADOS_WITH_QPDUNES
#include <acados/ocp_qp/ocp_qp_qpdunes.h>
#endif
#include <acados/dense_qp/dense_qp_hpipm.h>
#include <acados/dense_qp/dense_qp_qpoases.h>
#ifdef ACADOS_WITH_QORE
#include <acados/dense_qp/dense_qp_qore.h>
#endif

// #ifndef ACADOS_WITH_QPDUNES
#define ELIMINATE_X0
// #endif
#define GENERAL_CONSTRAINT_AT_TERMINAL_STAGE

#define NREP 100

#include "./mass_spring.c"

#define LAM_INIT 1.0

#define T_INIT 1.0

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" mass spring partial tightening example: acados + HPMPC\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
     * ocp qp
     ************************************************/

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    ocp_qp_dims *qp_dims = qp_in->dim;

    /************************************************
     * ocp qp solution
     ************************************************/

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_dims);

    /************************************************
     * ocp qp solvers
     ************************************************/

    // choose ocp qp solvers
    ocp_qp_solver_t ocp_qp_solvers[] =
    {
        PARTIAL_CONDENSING_HPMPC,
    };

    // choose values for N2 in partial condensing solvers
    int num_N2_values = 1;
    int num_M_values = 1;
    int N2_values[1]  = {15};
    int M_values[1]   = {10};

    int ii_max = 1;

    for (int ii = 0; ii < ii_max; ii++)
    {
        ocp_qp_solver_plan plan;
        plan.qp_solver = ocp_qp_solvers[ii];

        void *args = ocp_qp_create_args(&plan, qp_dims);

        for (int jj = 0; jj < num_N2_values; jj++) {
            for (int kk = 0; kk < num_M_values; kk++) {
                int N2 = N2_values[jj];
                int M = M_values[kk];

                printf("\nPartial condensing + partial tightening + HPMPC (N2 = %d, M = %d):\n\n", N2, M);
                ((ocp_qp_partial_condensing_args *)((ocp_qp_sparse_solver_args *)args)->pcond_args)->N2 = N2;
                ((ocp_qp_hpmpc_args *)((ocp_qp_sparse_solver_args *)args)->solver_args)->max_iter = 30;
                ((ocp_qp_hpmpc_args *)((ocp_qp_sparse_solver_args *)args)->solver_args)->M = M;
                
                ocp_qp_solver *qp_solver = ocp_qp_create(&plan, qp_dims, args);

                ocp_qp_hpmpc_memory *hpmpc_mem = (ocp_qp_hpmpc_memory *)((ocp_qp_sparse_solver_memory *)qp_solver->mem)->solver_memory;
                
                // initialize additional variables
                for (int ii = 0; ii <= 15; ii++) {
                    for (int jj = 0; jj < 2*(qp_dims->nb[ii] + qp_dims->ng[ii]); jj++) {
                        hpmpc_mem->lam0[ii].pa[jj] = LAM_INIT; 
                        hpmpc_mem->t0[ii].pa[jj] = T_INIT;
                    }
                }
                
                int acados_return = 0;

                ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
                ocp_qp_info min_info;

                // run QP solver NREP times and record min timings
                for (int rep = 0; rep < NREP; rep++)
                {
                    acados_return += ocp_qp_solve(qp_solver, qp_in, qp_out);

                    if (rep == 0)
                    {
                        min_info.num_iter = info->num_iter;
                        min_info.total_time = info->total_time;
                        min_info.condensing_time = info->condensing_time;
                        min_info.solve_QP_time = info->solve_QP_time;
                        min_info.interface_time = info->interface_time;
                    }
                    else
                    {
                        assert(min_info.num_iter == info->num_iter && "QP solver not cold started!");

                        if (info->total_time < min_info.total_time)
                            min_info.total_time = info->total_time;
                        if (info->condensing_time < min_info.condensing_time)
                            min_info.condensing_time = info->condensing_time;
                        if (info->solve_QP_time < min_info.solve_QP_time)
                            min_info.solve_QP_time = info->solve_QP_time;
                        if (info->interface_time < min_info.interface_time)
                            min_info.interface_time = info->interface_time;
                    }
                }

                /************************************************
                 * compute residuals
                 ************************************************/

                ocp_qp_res *qp_res = create_ocp_qp_res(qp_dims);
                ocp_qp_res_ws *res_ws = create_ocp_qp_res_ws(qp_dims);
                compute_ocp_qp_res(qp_in, qp_out, qp_res, res_ws);

                /************************************************
                 * print solution
                 ************************************************/

                 print_ocp_qp_out(qp_out);

                /************************************************
                 * print residuals
                 ************************************************/

    //			 print_ocp_qp_res(qp_res);

                /************************************************
                 * compute infinity norm of residuals
                 ************************************************/

                double res[4];
                compute_ocp_qp_res_nrm_inf(qp_res, res);
                double max_res = 0.0;
                for (int ii = 0; ii < 4; ii++)
                    max_res = (res[ii] > max_res) ? res[ii] : max_res;

                /************************************************
                 * print stats
                 ************************************************/

                printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

                print_ocp_qp_info(&min_info);

                /************************************************
                 * free memory
                 ************************************************/

                free(qp_solver);

                if (plan.qp_solver >= FULL_CONDENSING_HPIPM) break;
            }
        }

        free(args);
    }

    free(qp_in);
    free(qp_out);

    return 0;
}

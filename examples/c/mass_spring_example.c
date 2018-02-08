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

#ifndef ACADOS_WITH_QPDUNES
#define ELIMINATE_X0
#endif
#define GENERAL_CONSTRAINT_AT_TERMINAL_STAGE

#define NREP 100

#include "./mass_spring.c"

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" mass spring example: acados ocp_qp solvers\n");
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
        SPARSE_QP_HPIPM,
        #if ACADOS_WITH_HPMPC
        SPARSE_QP_HPMPC,
        #endif
        #if ACADOS_WITH_QPDUNES
        SPARSE_QP_QPDUNES,
        #endif
        PARTIAL_CONDENSING_HPIPM,
        #if ACADOS_WITH_HPMPC
        PARTIAL_CONDENSING_HPMPC,
        #endif
        #if ACADOS_WITH_QPDUNES
        PARTIAL_CONDENSING_QPDUNES,
        #endif
        FULL_CONDENSING_HPIPM,
        #ifdef ACADOS_WITH_QORE
        FULL_CONDENSING_QORE,
        #endif
        FULL_CONDENSING_QPOASES,
    };

    // choose values for N2 in partial condensing solvers
    int num_N2_values = 3;
    int N2_values[3] = {15, 5, 3};

    int ii_max = 9;
    #ifndef ACADOS_WITH_HPMPC
    ii_max--;
    #endif
    #ifndef ACADOS_WITH_QPDUNES
    ii_max--;
    #endif
    #ifndef ACADOS_WITH_QORE
    ii_max--;
    #endif

    for (int ii = 0; ii < ii_max; ii++)
    {
        ocp_qp_solver_config config;
        config.qp_solver = ocp_qp_solvers[ii];

        ocp_qp_solver_fcn_ptrs *fcn_ptrs = create_ocp_qp_solver_fcn_ptrs(&config, qp_dims);

        void *args = ocp_qp_create_args(fcn_ptrs, qp_dims);

        for (int jj = 0; jj < num_N2_values; jj++)
        {
            int N2 = N2_values[jj];

            // NOTE(nielsvd): needs to be implemented using the acados_c/options.h interface
#ifdef ACADOS_WITH_QPDUNES
            ocp_qp_qpdunes_args *qpdunes_args;
#endif
            switch (config.qp_solver)
            {
                case SPARSE_QP_HPIPM:
                    printf("\nHPIPM:\n\n");
                    ((ocp_qp_hpipm_args *) args)->hpipm_args->iter_max = 30;
                    break;
                case SPARSE_QP_HPMPC:
#ifdef ACADOS_WITH_HPMPC
                    printf("\nHPMPC:\n\n");
                    ((ocp_qp_hpmpc_args *)args)->max_iter = 30;
#endif
                    break;
                case SPARSE_QP_QPDUNES:
#ifdef ACADOS_WITH_QPDUNES
                    printf("\nQPDUNES:\n\n");
                    qpdunes_args = (ocp_qp_qpdunes_args *)args;
                    #ifdef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
                    qpdunes_args->stageQpSolver = QPDUNES_WITH_QPOASES;
                    #endif
                    qpdunes_args->warmstart = 0;
#endif
                    break;
                case PARTIAL_CONDENSING_HPIPM:
                    printf("\nPartial condensing + HPIPM (N2 = %d):\n\n", N2);
                    ((ocp_qp_partial_condensing_args *)((ocp_qp_sparse_solver_args *)args)->pcond_args)->N2 = N2;
                    ((ocp_qp_hpipm_args *)((ocp_qp_sparse_solver_args *)args)->solver_args)->hpipm_args->iter_max = 30;
                    break;
                case PARTIAL_CONDENSING_HPMPC:
#ifdef ACADOS_WITH_HPMPC
                    printf("\nPartial condensing + HPMPC (N2 = %d):\n\n", N2);
                    ((ocp_qp_partial_condensing_args *)((ocp_qp_sparse_solver_args *)args)->pcond_args)->N2 = N2;
                    ((ocp_qp_hpmpc_args *)((ocp_qp_sparse_solver_args *)args)->solver_args)->max_iter = 30;
#endif
                    break;
                case PARTIAL_CONDENSING_QPDUNES:
#ifdef ACADOS_WITH_QPDUNES
                    printf("\nPartial condensing + qpDUNES (N2 = %d):\n\n", N2);
                    #ifdef ELIMINATE_X0
                    assert(1==0 && "qpDUNES does not support ELIMINATE_X0 flag!");
                    #endif
                    ocp_qp_sparse_solver_args *solver_args = (ocp_qp_sparse_solver_args *)args;
                    qpdunes_args = (ocp_qp_qpdunes_args *)solver_args->solver_args;
                    #ifdef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
                    qpdunes_args->stageQpSolver = QPDUNES_WITH_QPOASES;
                    #endif
                    qpdunes_args->warmstart = 0;
                    ((ocp_qp_partial_condensing_args *)((ocp_qp_sparse_solver_args *)args)->pcond_args)->N2 = N2;
#endif
                    break;
                case FULL_CONDENSING_HPIPM:
                    printf("\nFull condensing + HPIPM:\n\n");
                    // default options
                    break;
                case FULL_CONDENSING_QORE:
#ifdef ACADOS_WITH_QORE
                    printf("\nFull condensing + QORE:\n\n");
                    // default options
                    break;
#endif
                case FULL_CONDENSING_QPOASES:
                    printf("\nFull condensing + QPOASES:\n\n");
                    // default options
                    break;
                case PARTIAL_CONDENSING_OOQP:
                    break;
            }

            ocp_qp_solver *qp_solver = ocp_qp_create(fcn_ptrs, qp_dims, args);

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

//			 print_ocp_qp_out(qp_out);

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

            // if (config.qp_solver >= FULL_CONDENSING_HPIPM) break;
        }
        free(args);
        free(fcn_ptrs);
    }

    free(qp_in);
    free(qp_out);

    return 0;
}

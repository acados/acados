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
#include <assert.h>
#include <stdlib.h>
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/ocp_qp/ocp_qp_sparse_solver.h"
#ifdef ACADOS_WITH_QPDUNES
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#else
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#endif
#include "acados/utils/create.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#ifndef ACADOS_WITH_QPDUNES
#define ELIMINATE_X0  // NOTE(dimitris): not supported by qpDUNES
#endif

#define NREP 1

#include "./mass_spring.c"


int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + partial condensing solvers\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
    * ocp qp
    ************************************************/

    // TODO(dimitris): write a print_ocp_qp function
    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;

    /************************************************
    * ocp qp solution
    ************************************************/

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_in->dim);

    /************************************************
    * partial condensing solvers
    ************************************************/

    #ifdef ACADOS_WITH_QPDUNES
    int num_of_solvers = 2;
    ocp_qp_solver_t pcond_solvers[2] = {QPDUNES, HPIPM};
    #else
    int num_of_solvers = 1;
    ocp_qp_solver_t pcond_solvers[1] = {HPIPM};
    #endif

    // different condensing arguments
    int num_of_cond_arguments = 3;
    int cond_arguments[3] = {3, 5, 15};

    for (int ii = 0; ii < num_of_solvers; ii++)
    {
        ocp_qp_solver_t solver = pcond_solvers[ii];
        ocp_qp_sparse_solver_args *arg = ocp_qp_sparse_solver_create_arguments(qp_in->dim, solver);

        for (int jj = 0; jj < num_of_cond_arguments; jj++)
        {
            // change partial condensing arguments
            arg->pcond_args->N2 = cond_arguments[jj];

            // set solver specific arguments and print solver name
            switch(solver)
            {
                #ifdef ACADOS_WITH_QPDUNES
                case QPDUNES:
                    printf("\nQPDUNES:\n\n");
                    ((ocp_qp_qpdunes_args *)(arg->solver_args))->stageQpSolver = QPDUNES_WITH_QPOASES;
                    break;
                #endif
                case HPIPM:
                    printf("\nHPIPM:\n\n");
                    ((ocp_qp_hpipm_args *)(arg->solver_args))->hpipm_args->iter_max = 21;
                    break;
                default:
                    printf("\nUnknown solver!\n\n");
                    exit(1);
            }

            // create partial condensing solver memory
            ocp_qp_sparse_solver_memory *mem =
            ocp_qp_sparse_solver_create_memory(qp_in->dim, arg);

            // create partial condensing solver workspace
            void *work = malloc(ocp_qp_sparse_solver_calculate_workspace_size(qp_in->dim, arg));

            int acados_return;  // 0 normal; 1 max iter

            acados_timer timer;
            acados_tic(&timer);

            ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;
            ocp_qp_info min_info;
            min_info.total_time = min_info.condensing_time = min_info.solve_QP_time = min_info.interface_time = 1e10;

            for (int rep = 0; rep < NREP; rep++) {
                acados_return = ocp_qp_sparse_solver(qp_in, qp_out, arg, mem, work);

                if (info->total_time < min_info.total_time) min_info.total_time = info->total_time;
                if (info->condensing_time < min_info.condensing_time) min_info.condensing_time = info->condensing_time;
                if (info->solve_QP_time < min_info.solve_QP_time) min_info.solve_QP_time = info->solve_QP_time;
                if (info->interface_time < min_info.interface_time) min_info.interface_time = info->interface_time;
                if (rep == 0)
                    min_info.num_iter = info->num_iter;
                #ifndef ACADOS_WITH_QPDUNES
                else
                    assert(min_info.num_iter == info->num_iter && "QP solver not cold started!");
                #endif
            }

            double time = acados_toc(&timer)/NREP;

            /************************************************
            * extract solution
            ************************************************/

            ocp_qp_dims *dims = qp_in->dim;

            colmaj_ocp_qp_out *sol;
            void *memsol = malloc(colmaj_ocp_qp_out_calculate_size(dims));
            assign_colmaj_ocp_qp_out(dims, &sol, memsol);
            convert_ocp_qp_out_to_colmaj(qp_out, sol);

            /************************************************
            * compute residuals
            ************************************************/

            ocp_qp_res *qp_res = create_ocp_qp_res(dims);
            ocp_qp_res_ws *res_ws = create_ocp_qp_res_ws(dims);
            compute_ocp_qp_res(qp_in, qp_out, qp_res, res_ws);

            /************************************************
            * compute infinity norm of residuals
            ************************************************/

            double res[4];
            compute_ocp_qp_res_nrm_inf(qp_res, res);
            double max_res = 0.0;
            for (int ii = 0; ii < 4; ii++) max_res = (res[ii] > max_res) ? res[ii] : max_res;

            // solver specific checks
            if (solver == HPIPM) assert(max_res <= 1e6*ACADOS_EPS && "The largest KKT residual greater than 1e6*ACADOS_EPS");
            // else assert(max_res <= 1e1*ACADOS_EPS && "The largest KKT residual greater than ACADOS_EPS");

            /************************************************
            * print stats
            ************************************************/

            printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

            printf("\nN2 = %d\n\n\n", arg->pcond_args->N2);

            print_ocp_qp_info(&min_info);

            /************************************************
            * free memory
            ************************************************/

            // ocp_qp_qpdunes_free_memory(mem);
            free(sol);
            free(qp_res);
            free(res_ws);
            free(mem);
            free(work);
        }
        free(arg);
    }

    free(qp_in);
    free(qp_out);

    return 0;
}

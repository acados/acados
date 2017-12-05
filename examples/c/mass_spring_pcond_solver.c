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

// #define ELIMINATE_X0
#define NREP 100

#include "./mass_spring.c"


int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + partial condensing solver\n");
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
    * partial condensing solver
    ************************************************/

    // choose QP solver
    ocp_qp_solver_t qp_solver_name = HPIPM;

    // create partial condensing solver args
    #ifdef ACADOS_WITH_QPDUNES
    ocp_qp_sparse_solver_args *arg = ocp_qp_sparse_solver_create_arguments(qp_in->dim, QPDUNES);
    #else
    ocp_qp_sparse_solver_args *arg = ocp_qp_sparse_solver_create_arguments(qp_in->dim, HPIPM);
    #endif

    // change partial condensing arguments
    arg->pcond_args->N2 = 5;

    // change qp solver arguments
    #ifdef ACADOS_WITH_QPDUNES
    ((ocp_qp_qpdunes_args *)(arg->solver_args))->stageQpSolver = QPDUNES_WITH_QPOASES;
    #else
    ((ocp_qp_hpipm_args *)(arg->solver_args))->hpipm_args->iter_max = 21;
    printf("maxIter = %d\n", ((ocp_qp_hpipm_args *)(arg->solver_args))->hpipm_args->iter_max);
    #endif

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
    * print solution and stats
    ************************************************/

    printf("\nu = \n");
    for (int ii = 0; ii < N; ii++) d_print_mat(1, nu[ii], sol->u[ii], 1);

    printf("\nx = \n");
    for (int ii = 0; ii <= N; ii++) d_print_mat(1, nx[ii], sol->x[ii], 1);

    printf("\npi = \n");
    for (int ii = 0; ii < N; ii++) d_print_mat(1, nx[ii+1], sol->pi[ii], 1);

    printf("\nlam = \n");
    for (int ii = 0; ii <= N; ii++) d_print_mat(1, 2*nb[ii]+2*ng[ii], sol->lam[ii], 1);

    printf("\nSolution time for N2 = %d, averaged over %d runs: %5.2e seconds\n\n\n",
        arg->pcond_args->N2, NREP, time);

    print_ocp_qp_info(&min_info);

    /************************************************
    * free memory
    ************************************************/

    free(qp_in);
    free(qp_out);
    free(sol);
    free(arg);
    free(mem);
    free(work);

    return 0;
}

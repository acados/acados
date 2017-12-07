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
#include <acados/ocp_qp/ocp_qp_sparse_solver.h>
#include <acados/ocp_qp/ocp_qp_full_condensing_solver.h>
#include <acados/ocp_qp/ocp_qp_hpipm.h>
#include <acados/dense_qp/dense_qp_hpipm.h>
#include <acados/dense_qp/dense_qp_qpoases.h>
#include <acados/dense_qp/dense_qp_qore.h>

#define ELIMINATE_X0
#define NREP 100

#include "./mass_spring.c"

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + hpipm\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
     * ocp qp
     ************************************************/

    // TODO(dimitris): write a print_ocp_qp function
    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    ocp_qp_dims *qp_dims = qp_in->dim;

    int N = qp_dims->N;
    int *nx = qp_dims->nx;
    int *nu = qp_dims->nu;
    int *nb = qp_dims->nb;
    int *ng = qp_dims->ng;

    /************************************************
     * ocp qp solution
     ************************************************/

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_dims);

    /************************************************
     * ipm
     ************************************************/
    ocp_qp_solver_plan plan;
    plan.qp_solver = FULL_CONDENSING_QPOASES;

    void *args = ocp_qp_create_args(&plan, qp_dims);

    // NOTE(nielsvd): needs to be implemented using the acados_c/options.h interface
    switch (plan.qp_solver) {
        case PARTIAL_CONDENSING_HPIPM:
            ((ocp_qp_partial_condensing_args *)((ocp_qp_sparse_solver_args *)args)->pcond_args)->N2 = N;
            ((ocp_qp_hpipm_args *)((ocp_qp_sparse_solver_args *)args)->solver_args)->hpipm_args->iter_max = 10;
            break;
        case FULL_CONDENSING_HPIPM:
            // default options
            break;
        case FULL_CONDENSING_QORE:
            // default options
            break;
        case FULL_CONDENSING_QPOASES:
            // default options
            break;
    }

    ocp_qp_solver *qp_solver = ocp_qp_create(&plan, qp_dims, args);

    int acados_return;  // 0 normal; 1 max iter

    acados_timer timer;
    acados_tic(&timer);

    for (int rep = 0; rep < NREP; rep++) {
        acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);
    }

    double time = acados_toc(&timer) / NREP;

    /************************************************
     * extract solution
     ************************************************/

    colmaj_ocp_qp_out *sol;
    void *memsol = malloc(colmaj_ocp_qp_out_calculate_size(qp_dims));
    assign_colmaj_ocp_qp_out(qp_dims, &sol, memsol);
    convert_ocp_qp_out_to_colmaj(qp_out, sol);

    /************************************************
     * compute residuals
     ************************************************/

    ocp_qp_res *qp_res = create_ocp_qp_res(qp_dims);
    ocp_qp_res_ws *res_ws = create_ocp_qp_res_ws(qp_dims);
    compute_ocp_qp_res(qp_in, qp_out, qp_res, res_ws);

    /************************************************
     * compute infinity norm of residuals
     ************************************************/

    double res[4];
    compute_ocp_qp_res_nrm_inf(qp_res, res);
    double max_res = 0.0;
    for (int ii = 0; ii < 4; ii++)
        max_res = (res[ii] > max_res) ? res[ii] : max_res;
    assert(max_res <= 1e6 * ACADOS_EPS &&
           "The largest KKT residual greater than 1e6*ACADOS_EPS");

    /************************************************
     * print solution and stats
     ************************************************/

    printf("\nu = \n");
    for (int ii = 0; ii < N; ii++) d_print_mat(1, nu[ii], sol->u[ii], 1);

    printf("\nx = \n");
    for (int ii = 0; ii <= N; ii++) d_print_mat(1, nx[ii], sol->x[ii], 1);

    printf("\npi = \n");
    for (int ii = 0; ii < N; ii++) d_print_mat(1, nx[ii + 1], sol->pi[ii], 1);

    printf("\nlam = \n");
    for (int ii = 0; ii <= N; ii++)
        d_print_mat(1, 2 * nb[ii] + 2 * ng[ii], sol->lam[ii], 1);

    // NOTE(nielsvd): how can we improve/generalize this?
    ocp_qp_hpipm_memory *pcond_hpipm_mem = (ocp_qp_hpipm_memory *)((ocp_qp_sparse_solver_memory *) qp_solver->mem)->solver_memory;
    dense_qp_hpipm_memory *fcond_hpipm_mem = (dense_qp_hpipm_memory *)((ocp_qp_full_condensing_solver_memory *) qp_solver->mem)->solver_memory;
    dense_qp_qpoases_memory *fcond_qpoases_mem = (dense_qp_qpoases_memory *)((ocp_qp_full_condensing_solver_memory *) qp_solver->mem)->solver_memory;
    switch (plan.qp_solver) {
        case PARTIAL_CONDENSING_HPIPM:
            printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2],
                   res[3]);
            printf("\nSolution time for %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n",
                    pcond_hpipm_mem->hpipm_workspace->iter, NREP, time);
            break;
        case FULL_CONDENSING_HPIPM:
            printf("\nSolution time for %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n", fcond_hpipm_mem->hpipm_workspace->iter, NREP, time);
            break;
        case FULL_CONDENSING_QORE:
            // no post-processing
            break;
        case FULL_CONDENSING_QPOASES:
            printf("\nSolution time for %d NWSR, averaged over %d runs: %5.2e seconds\n\n\n", fcond_qpoases_mem->nwsr, NREP, time);
            break;
    }

    /************************************************
     * free memory
     ************************************************/

    free(qp_in);
    free(qp_out);
    free(sol);
    free(qp_solver);
    free(args);

    return 0;
}

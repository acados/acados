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
// acados_c
#include <acados_c/ocp_qp.h>
#include <acados_c/legacy_create.h>
// acados
#include <acados/ocp_qp/ocp_qp_partial_condensing.h>
// NOTE(nielsvd): required to cast memory etc. should go.
#include <acados/ocp_qp/ocp_qp_sparse_solver.h>
#include <acados/ocp_qp/ocp_qp_hpipm.h>

#include <acados/utils/print.h>

#define ELIMINATE_X0
#define NREP 100

#include "./mass_spring.c"


int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + partial condensing + hpipm\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
    * ocp qp
    ************************************************/

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    ocp_qp_dims *qp_dims = qp_in->dim;

    int N = qp_dims->N;
    int *nx = qp_dims->nx;
    int *nu = qp_dims->nu;
    int *nb = qp_dims->nb;
    int *ng = qp_dims->ng;

    /************************************************
    * partial condensing arguments/memory
    ************************************************/


    ocp_qp_partial_condensing_args *pcond_args =
        ocp_qp_partial_condensing_create_arguments(qp_in->dim, NULL);

    pcond_args->N2 = 4;

    ocp_qp_partial_condensing_memory *pcond_mem =
        ocp_qp_partial_condensing_create_memory(qp_in->dim, pcond_args);

    /************************************************
    * partially condensed ocp qp
    ************************************************/

    ocp_qp_in *pcond_qp_in = create_ocp_qp_in(pcond_args->pcond_dims);

    // print_ocp_qp_dims(qp_in->dim);
    // print_ocp_qp_dims(pcond_qp_in->dim);

    /************************************************
    * ocp qp solution
    ************************************************/

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_dims);

    /************************************************
    * partially condensed ocp qp solution
    ************************************************/

    ocp_qp_out *pcond_qp_out = create_ocp_qp_out(pcond_qp_in->dim);

    /************************************************
    * ipm
    ************************************************/

    ocp_qp_solver_config config;
    config.qp_solver = PARTIAL_CONDENSING_HPIPM;

    ocp_qp_solver_fcn_ptrs *fcn_ptrs = create_ocp_qp_solver_fcn_ptrs(&config, qp_dims);

    void *arg = ocp_qp_create_args(fcn_ptrs, qp_dims);

    // NOTE(nielsvd): needs to be implemented using the acados_c/options.h interface
    ((ocp_qp_partial_condensing_args *)((ocp_qp_sparse_solver_args *)arg)->pcond_args)->N2 = pcond_args->N2;
    ((ocp_qp_hpipm_args *)((ocp_qp_sparse_solver_args *)arg)->solver_args)->hpipm_args->iter_max = 10;
    ((ocp_qp_hpipm_args *)((ocp_qp_sparse_solver_args *)arg)->solver_args)->hpipm_args->res_g_max = 1e-8;

    ocp_qp_solver *qp_solver = ocp_qp_create(fcn_ptrs,pcond_args->pcond_dims, arg);

    free(fcn_ptrs);

	int acados_return;  // 0 normal; 1 max iter

    acados_timer timer;
    acados_tic(&timer);

	for (int rep = 0; rep < NREP; rep++) {

        ocp_qp_partial_condensing(qp_in, pcond_qp_in, pcond_args, pcond_mem, NULL);

        acados_return = ocp_qp_solve(qp_solver, pcond_qp_in, pcond_qp_out);

        ocp_qp_partial_expansion(pcond_qp_out, qp_out, pcond_args, pcond_mem, NULL);
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
    assert(max_res <= 1e6*ACADOS_EPS && "The largest KKT residual greater than 1e6*ACADOS_EPS");

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

    // NOTE(nielsvd): how can we improve/generalize this?
    ocp_qp_hpipm_memory *mem = (ocp_qp_hpipm_memory *)((ocp_qp_sparse_solver_memory *) qp_solver->mem)->solver_memory;

    printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

    printf("\nNumber of %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        mem->hpipm_workspace->iter, NREP, time);

    /************************************************
    * free memory
    ************************************************/

    free(qp_in);
    free(pcond_qp_in);
    free(qp_out);
    free(pcond_qp_out);
    free(sol);
    free(qp_solver);
    free(qp_res);
    free(res_ws);
    free(arg);
    free(pcond_args);
    free(pcond_mem);
    return 0;
}

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
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/ocp_qp/ocp_qp_condensing_solver.h"
#include "acados/utils/create.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#define ELIMINATE_X0
#define NREP 1000

#include "./mass_spring.c"

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + condensing_hpipm\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
    * ocp qp
    ************************************************/

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *ns = qp_in->dim->ns;

    /************************************************
    * ocp qp solution
    ************************************************/

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_in->dim);

    /************************************************
    * ipm
    ************************************************/

    ocp_qp_condensing_solver_args *arg = ocp_qp_condensing_solver_create_arguments(qp_in->dim, CONDENSING_HPIPM);

    ocp_qp_condensing_solver_memory *mem = ocp_qp_condensing_solver_create_memory(qp_in->dim, arg);

    void *work = malloc(ocp_qp_condensing_solver_calculate_workspace_size(qp_in->dim, arg));

	int acados_return; // 0 normal; 1 max iter

    acados_timer timer;
    acados_tic(&timer);

	for (int rep = 0; rep < NREP; rep++) {
        acados_return = ocp_qp_condensing_solver(qp_in, qp_out, arg, mem, work);
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

    ocp_qp_res *qp_res = create_ocp_qp_res(qp_in->dim);
    ocp_qp_res_ws *res_ws = create_ocp_qp_res_ws(qp_in->dim);
    compute_ocp_qp_res(qp_in, qp_out, qp_res, res_ws);

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

    printf("\nres_g = \n");
    for (int ii = 0; ii <= N; ii++) d_print_tran_strvec(nu[ii]+nx[ii]+2*ns[ii], qp_res->res_g+ii, 0);

    printf("\nres_b = \n");
    for (int ii = 0; ii < N; ii++) d_print_tran_strvec(nx[ii+1], qp_res->res_b+ii, 0);

    printf("\nres_d = \n");
    for (int ii = 0; ii <= N; ii++) d_print_tran_strvec(2*nb[ii]+2*ng[ii]+2*ns[ii], qp_res->res_d+ii, 0);

    printf("\nres_m = \n");
    for (int ii = 0; ii <= N; ii++) d_print_tran_strvec(2*nb[ii]+2*ng[ii]+2*ns[ii], qp_res->res_m+ii, 0);


    dense_qp_hpipm_memory *tmp_mem = (dense_qp_hpipm_memory *) mem->solver_memory;

    printf("\nSolution time for %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        tmp_mem->hpipm_workspace->iter, NREP, time);

    /************************************************
    * free memory
    ************************************************/

    free(qp_in);
    free(qp_out);
    free(sol);
    free(arg);
    free(mem);

	return 0;
}

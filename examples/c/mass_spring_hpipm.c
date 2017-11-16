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
#include "acados/ocp_qp/ocp_qp_hpipm.h"
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
    printf(" acados + hpipm\n");
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
    * ipm
    ************************************************/

    ocp_qp_hpipm_args *arg = ocp_qp_hpipm_create_arguments(qp_in->dim);

    arg->hpipm_args->iter_max = 10;

    ocp_qp_hpipm_memory *mem = ocp_qp_hpipm_create_memory(qp_in->dim, arg);

	int acados_return;  // 0 normal; 1 max iter

    acados_timer timer;
    acados_tic(&timer);

	for (int rep = 0; rep < NREP; rep++) {
        acados_return = ocp_qp_hpipm(qp_in, qp_out, arg, mem);
	}

    double time = acados_toc(&timer)/NREP;

    /************************************************
    * extract solution
    ************************************************/

    // ocp_qp_dims dims;
    // dims.N = N;
    // dims.nx = nx;
    // dims.nu = nu;
    // dims.nb = nb;
    // dims.ns = ns;
    // dims.ng = ng;

    ocp_qp_dims *dims = qp_in->dim;

    col_maj_ocp_qp_out *sol;
    void *memsol = malloc(col_maj_ocp_qp_out_calculate_size(dims));
    assign_col_maj_ocp_qp_out(dims, &sol, memsol);
    convert_to_col_maj_ocp_qp_out(dims, qp_out, sol);

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

    printf("\ninf norm res: %e, %e, %e, %e, %e\n", mem->hpipm_workspace->qp_res[0],
           mem->hpipm_workspace->qp_res[1], mem->hpipm_workspace->qp_res[2],
           mem->hpipm_workspace->qp_res[3], mem->hpipm_workspace->res_workspace->res_mu);

    printf("\nSolution time for %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        mem->hpipm_workspace->iter, NREP, time);

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

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
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/utils/create.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#define NREP 1000
#define ELIMINATE_X0

#include "./mass_spring.c"

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + condensing + hpipm + expansion\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
    * ocp qp
    ************************************************/

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    int N = qp_in->size->N;
    int *nx = qp_in->size->nx;
    int *nu = qp_in->size->nu;
    int *nb = qp_in->size->nb;
    int *ng = qp_in->size->ng;

    /************************************************
    * dense qp
    ************************************************/

    // dummy dense qp to calculate dimensions
    dense_qp_in qpd_dummy;
    dummy_dense_qp_in(&qpd_dummy, qp_in->size);

    int nvd = qpd_dummy.nv;
    int ned = qpd_dummy.ne;
    int ngd = qpd_dummy.ng;
    int nbd = qpd_dummy.nb;
    int nsd = qpd_dummy.ns;

    dense_qp_in *qpd_in = create_dense_qp_in(nvd, ned, nbd, ngd, nsd);

    ocp_qp_condensing_args *cond_args = ocp_qp_condensing_create_arguments(qp_in->size);
    ocp_qp_condensing_memory *cond_memory = ocp_qp_condensing_create_memory(qp_in->size, cond_args);

    /************************************************
    * ocp qp solution
    ************************************************/

    ocp_qp_out *qp_out = create_ocp_qp_out(qp_in->size);

    /************************************************
    * dense sol
    ************************************************/

    dense_qp_out *qpd_out = create_dense_qp_out(nvd, ned, nbd, ngd, nsd);

    /************************************************
    * dense ipm
    ************************************************/

    dense_qp_hpipm_args *argd = dense_qp_hpipm_create_arguments(qpd_in);

    // argd->hpipm_args->iter_max = 10;

    dense_qp_hpipm_memory *mem = dense_qp_hpipm_create_memory(qpd_in, argd);
	int acados_return;  // 0 normal; 1 max iter

    acados_timer timer;
    acados_tic(&timer);

	for(int rep = 0; rep < NREP; rep++) {

        ocp_qp_condensing(qp_in, qpd_in, cond_args, cond_memory);

        acados_return = dense_qp_hpipm(qpd_in, qpd_out, argd, mem);

        ocp_qp_expansion(qpd_out, qp_out, cond_args, cond_memory);
    }

    real_t time = acados_toc(&timer)/NREP;

    /************************************************
    * extract solution
    ************************************************/

    ocp_qp_dims *dims = qp_in->size;

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

    printf("\nSolution time for %d IPM iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        mem->hpipm_workspace->iter, NREP, time);

    /************************************************
    * free memory
    ************************************************/

    free(qp_in);
    free(qpd_in);
    free(qp_out);
    free(qpd_out);
    free(sol);
    free(argd);
    free(mem);
    free(cond_memory);

	return 0;
}



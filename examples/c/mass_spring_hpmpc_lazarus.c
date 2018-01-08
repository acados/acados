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

 #define _GNU_SOURCE
 #include <fenv.h>
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
#define NREP 100

#include "./mass_spring.c"


int main() {
    // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + hpmpc\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
    * ocp qp
    ************************************************/

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring();

    int N   = qp_in->dim->N;
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

    ocp_qp_hpmpc_args *arg = ocp_qp_hpmpc_create_arguments(qp_in->dim);

    // arg->hpmpc_args->iter_max = 10;

    ocp_qp_hpmpc_memory *mem = ocp_qp_hpmpc_create_memory(qp_in->dim, qp_in, arg);

    // intialize variables
    for (int ii = 0; ii <= N;  ii++) {
        for (int jj = 0; jj < 2*nb[ii] + 2*ng[ii]; jj++) {
            arg->t0[ii][jj] = 10;
            arg->lam0[ii][jj] = 10;
        }

    }
	int acados_return;  // 0 normal; 1 max iter

    acados_timer timer;
    acados_tic(&timer);

	for (int rep = 0; rep < NREP; rep++) {
        acados_return = ocp_qp_hpmpc(qp_in, qp_out, arg, mem, NULL);
	}

    double time = acados_toc(&timer)/NREP;

    printf("acados_return = %i", acados_return);

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

    // printf("\ninf norm res: %e, %e, %e, %e, %e\n", mem->hpmpc_workspace->qp_res[0],
        //    mem->hpmpc_workspace->qp_res[1], mem->hpmpc_workspace->qp_res[2],
        //    mem->hpmpc_workspace->qp_res[3], mem->hpmpc_workspace->res_workspace->res_mu);

    printf("\nN ITER = %i\n\n\n",
        arg->out_iter);

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

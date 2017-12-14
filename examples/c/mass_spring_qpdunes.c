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

#ifdef ACADOS_WITH_QPDUNES

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/utils/create.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"

#define NREP 100

#include "./mass_spring.c"


int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" acados + qpDUNES\n");
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
    * qpDUNES
    ************************************************/

    ocp_qp_qpdunes_args *arg = ocp_qp_qpdunes_create_arguments(qp_in->dim);

    arg->stageQpSolver = QPDUNES_WITH_QPOASES;

    ocp_qp_qpdunes_memory *mem = ocp_qp_qpdunes_create_memory(qp_in->dim, arg);

    int size_work = ocp_qp_qpdunes_calculate_workspace_size(qp_in->dim, arg);
    void *work = acados_malloc(size_work, 1);;

	int acados_return;

    acados_timer timer;
    acados_tic(&timer);

	for (int rep = 0; rep < NREP; rep++) {
        acados_return = ocp_qp_qpdunes(qp_in, qp_out, arg, mem, work);
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

    printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

    printf("\nSolution time for %d dual Newton iterations, averaged over %d runs: %5.2e seconds\n\n\n",
        -1, NREP, time);

    /************************************************
    * free memory
    ************************************************/

    // TODO(dimitris): there are still memory leaks from qpOASES
    ocp_qp_qpdunes_free_memory(mem);
    free(qp_in);
    free(qp_out);
    free(qp_res);
    free(res_ws);
    free(arg);
    free(mem);
    free(work);
    free(sol);

    return 0;
}

#else

int main( )
{
    printf("qpDUNES example skipped\n");
}

#endif
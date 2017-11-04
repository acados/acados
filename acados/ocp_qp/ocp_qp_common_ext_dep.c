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
#include <stdlib.h>
#include <assert.h>
#include <string.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
// acados
#include "acados/ocp_qp/ocp_qp_common_ext_dep.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

// TODO TEMP
#include "acados/ocp_qp/ocp_qp_hpipm_ext_dep.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases_ext_dep.h"


ocp_qp_in *create_ocp_qp_in(ocp_qp_dims *dims) {

    ocp_qp_in *qp_in;

    int size = ocp_qp_in_calculate_size(dims);
    void *ptr = malloc(size);
    char *ptr_end = assign_ocp_qp_in(dims, &qp_in, ptr);
    assert((char*)ptr + size >= ptr_end); (void) ptr_end;

    return qp_in;
}



ocp_qp_out *create_ocp_qp_out(ocp_qp_dims *dims) {

    ocp_qp_out *qp_out;

    int size = ocp_qp_out_calculate_size(dims);
    void *ptr = malloc(size);
    char *ptr_end = assign_ocp_qp_out(dims, &qp_out, ptr);
    assert((char*)ptr + size >= ptr_end); (void) ptr_end;

    return qp_out;
}



void print_ocp_qp_dims(ocp_qp_dims *dims) {

    int N = dims->N;

    printf("k\tnx\tnu\tnb\tnbx\tnbu\tng\tns\n");

    for (int kk = 0; kk < N+1; kk++) {
        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n", kk, dims->nx[kk], dims->nu[kk], dims->nb[kk],
            dims->nbx[kk], dims->nbu[kk], dims->ng[kk], dims->ns[kk]);
    }

    printf("\nmemsize = %d\n", dims->memsize);

}


void print_ocp_qp_in(ocp_qp_in *qp_in)
{
    int N = qp_in->size->N;
    int *nx = qp_in->size->nx;
    int *nu = qp_in->size->nu;
    int *nb = qp_in->size->nb;
    int *ng = qp_in->size->ng;

    for (int ii = 0; ii < N+1; ii++)
    {
        printf("k = %d\n\n", ii);

        printf("RSQrq =\n");
        d_print_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &qp_in->RSQrq[ii], 0 , 0);

        printf("rq =\n");
        d_print_tran_strvec(nu[ii]+nx[ii], &qp_in->rq[ii], 0);


        if (ii < N)
        {
            printf("BAbt =\n");
            d_print_strmat(nu[ii]+nx[ii]+1, nx[ii+1], &qp_in->BAbt[ii], 0 , 0);

            printf("b =\n");
            d_print_tran_strvec(nx[ii+1], &qp_in->b[ii], 0);
        }

        printf("idxb = %d\n", qp_in->idxb[ii][0]);
        int_print_mat(1, nb[ii], qp_in->idxb[ii], 1);

        printf("d =\n");
        d_print_tran_strvec(2*nb[ii]+2*ng[ii], &qp_in->d[ii], 0);
    }

}
// TEMP TO FIX GN_SQP

// void ocp_qp_in_copy_dynamics(const real_t *A, const real_t *B, const real_t *b,
//     ocp_qp_in *qp_in, int_t stage) {

//     real_t **hA = (real_t **) qp_in->A;
//     real_t **hB = (real_t **) qp_in->B;
//     real_t **hb = (real_t **) qp_in->b;

//     memcpy(hA[stage], A, qp_in->nx[stage+1]*qp_in->nx[stage]*sizeof(real_t));
//     memcpy(hB[stage], B, qp_in->nx[stage+1]*qp_in->nu[stage]*sizeof(real_t));
//     memcpy(hb[stage], b, qp_in->nx[stage+1]*sizeof(real_t));
// }

// void ocp_qp_in_copy_objective(const real_t *Q, const real_t *S, const real_t *R, const real_t *q,
//      const real_t *r, ocp_qp_in *qp_in, int_t stage) {

//     real_t **hQ = (real_t **) qp_in->Q;
//     real_t **hS = (real_t **) qp_in->S;
//     real_t **hR = (real_t **) qp_in->R;
//     real_t **hq = (real_t **) qp_in->q;
//     real_t **hr = (real_t **) qp_in->r;

//     memcpy(hQ[stage], Q, qp_in->nx[stage]*qp_in->nx[stage]*sizeof(real_t));
//     memcpy(hS[stage], S, qp_in->nu[stage]*qp_in->nx[stage]*sizeof(real_t));
//     memcpy(hR[stage], R, qp_in->nu[stage]*qp_in->nu[stage]*sizeof(real_t));
//     memcpy(hq[stage], q, qp_in->nx[stage]*sizeof(real_t));
//     memcpy(hr[stage], r, qp_in->nu[stage]*sizeof(real_t));
// }

ocp_qp_solver *create_ocp_qp_solver(const ocp_qp_in *qp_in, const char *solver_name,
    void *solver_options) {
    ocp_qp_solver *qp_solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver));

    qp_solver->qp_in = (ocp_qp_in *) qp_in;
    qp_solver->qp_out = create_ocp_qp_out(qp_in->size);
    qp_solver->args = solver_options;

    if (!strcmp(solver_name, "qpdunes")) {
    // if (qp_solver->args == NULL)
    //     qp_solver->args = ocp_qp_qpdunes_create_arguments(QPDUNES_NONLINEAR_MPC);
    // qp_solver->fun = &ocp_qp_qpdunes;
    // qp_solver->initialize = &ocp_qp_qpdunes_initialize;
    // qp_solver->destroy = &ocp_qp_qpdunes_destroy;
    #ifdef OOQP
    // } else if (!strcmp(solver_name, "ooqp")) {
    // if (qp_solver->args == NULL)
    //     qp_solver->args = ocp_qp_ooqp_create_arguments();
    // qp_solver->fun = &ocp_qp_ooqp;
    // qp_solver->initialize = &ocp_qp_ooqp_initialize;
    // qp_solver->destroy = &ocp_qp_ooqp_destroy;
    #endif
    } else if (!strcmp(solver_name, "condensing_qpoases")) {
    if (qp_solver->args == NULL)
        qp_solver->args = ocp_qp_condensing_qpoases_create_arguments(qp_in->size);
    qp_solver->fun = &ocp_qp_condensing_qpoases;
    qp_solver->initialize = &ocp_qp_condensing_qpoases_initialize;
    qp_solver->destroy = &ocp_qp_condensing_qpoases_destroy;
    } else if (!strcmp(solver_name, "condensing_hpipm")) {
    // if (qp_solver->args == NULL)
    //     qp_solver->args = ocp_qp_condensing_hpipm_create_arguments(qp_in);
    // qp_solver->fun = &ocp_qp_condensing_hpipm;
    // qp_solver->initialize = &ocp_qp_condensing_hpipm_initialize;
    // qp_solver->destroy = &ocp_qp_condensing_hpipm_destroy;
    } else if (!strcmp(solver_name, "hpipm")) {
    if (qp_solver->args == NULL)
        qp_solver->args = ocp_qp_hpipm_create_arguments(qp_in->size);
    qp_solver->fun = &ocp_qp_hpipm;
    qp_solver->initialize = &ocp_qp_hpipm_initialize;
    qp_solver->destroy = &ocp_qp_hpipm_destroy;
    } else {
    printf("Chosen QP solver not available\n");
    exit(1);
    }
    qp_solver->initialize(qp_solver->qp_in, qp_solver->args, &qp_solver->mem, &qp_solver->work);

    return qp_solver;
}
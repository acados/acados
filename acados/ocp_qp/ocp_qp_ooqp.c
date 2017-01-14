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

#include <stdio.h>
#include <stdlib.h>
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "acados/utils/print.h"
#include "OOQP/include/cQpGenSparse.h"
// #include "blasfeo/include/blasfeo_i_aux.h"

int_t ocp_qp_ooqp_create_workspace(const ocp_qp_in *in, ocp_qp_ooqp_workspace *work) {
    int_t ii;
    int N = in->N;
    work->nx = 0;    // # of primal optimization variables
    work->nnzQ = 0;  // # non-zeros in lower part of Hessian
    work->nnzA = 0;  // TODO(dimitris): # non-zeros in matrix of equality constraints
    work->nnzC = 0;  // TODO(dimitris): # non-zeros in matrix of inequality constraints
    work->my = 0;    // # of equality constraints
    work->mz = 0;    // # of inequality constraints

    for (ii = 0; ii < N; ii++) {
        work->nx += in->nx[ii] + in->nu[ii];
        work->nnzQ += (in->nx[ii]*in->nx[ii] - in->nx[ii])/2 + in->nx[ii];
        work->nnzQ += (in->nu[ii]*in->nu[ii] - in->nu[ii])/2 + in->nu[ii];
        work->nnzQ += in->nx[ii]*in->nu[ii];
        work->my += N*in->nx[ii+1]; // TODO(dimitris): double check this
        work->mz += in->nb[ii];
    }
    work->nx += in->nx[N];
    work->nnzQ += (in->nx[N]*in->nx[N] -in->nx[N])/2 + in->nx[N];
    work->mz += in->nb[N];

    // TODO(dimitris): Check with giaf how x0 constr. is handled

    // initialization
    printf("nnzQ = %d\n", work->nnzQ);
    newQpGenSparse(&work->c, work->nx,
        &work->irowQ, work->nnzQ, &work->jcolQ, &work->dQ,
        &work->xlow, &work->ixlow, &work->xupp, &work->ixupp,
        &work->irowA, work->nnzA, &work->jcolA, &work->dA,
        &work-> bA, work->my,
        &work->irowC, work->nnzC, &work->jcolC, &work->dC,
        &work->clow, work->mz, &work->iclow, &work->cupp, &work->icupp, &work->ierr);
    // printf("flag = %d\n", work->ierr);
}

int_t ocp_qp_ooqp(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
    void *args_, void *workspace_) {

int_t return_flag;

return return_flag;
}

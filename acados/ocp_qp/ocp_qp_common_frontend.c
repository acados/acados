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
// TODO(dimitris): remove memcpy to avoid this dependency?
#include <string.h>
#include <assert.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
// hpipm
// #include "hpipm_d_ocp_qp.h"
// #include "hpipm_d_ocp_qp_sol.h"
// acados
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"
// #include "utils/types.h"
// #include "ocp_qp/ocp_qp_hpipm.h"


int col_maj_ocp_qp_in_calculate_size(ocp_qp_dims *dims)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *nc = dims->ng;

    int size = sizeof(col_maj_ocp_qp_in);

    size += 4*(N+1)*sizeof(int);  // nx, nu, nb, nc
    size += 3*N*sizeof(double *);  // A, B, b
    size += 11*(N+1)*sizeof(double *);  // ...
    size += 1*(N+1)*sizeof(int *);  // idxb

    for (int k = 0; k < N+1; k++) {

    if (k < N) {
    size += nx[k+1]*nx[k]*sizeof(double);  // A
    size += nx[k+1]*nu[k]*sizeof(double);  // B
    size += nx[k+1]*sizeof(double);  // b
    }

    size += nx[k]*nx[k]*sizeof(double);  // Q
    size += nu[k]*nx[k]*sizeof(double);  // S
    size += nu[k]*nu[k]*sizeof(double);  // R
    size += nx[k]*sizeof(double);  // q
    size += nu[k]*sizeof(double);  // r
    size += nb[k]*sizeof(int);  // idxb
    size += 2*nb[k]*sizeof(double);  // lb, ub
    size += nc[k]*nx[k]*sizeof(double);  // Cx
    size += nc[k]*nu[k]*sizeof(double);  // Cu
    size += 2*nc[k]*sizeof(double);  // lc, uc
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *col_maj_ocp_qp_in_assign(ocp_qp_dims *dims, col_maj_ocp_qp_in **qp_in, void *ptr)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *nc = dims->ng;

    // pointer to initialize QP data to zero
    char *c_ptr_QPdata;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_in = (col_maj_ocp_qp_in *) c_ptr;
    c_ptr += sizeof(col_maj_ocp_qp_in);

    // copy dimensions to workspace
    (*qp_in)->N = N;

    (*qp_in)->nx = (int *) c_ptr;
    memcpy(c_ptr, nx, (N+1)*sizeof(int));
    c_ptr += (N+1)*sizeof(int);

    (*qp_in)->nu = (int *) c_ptr;
    memcpy(c_ptr, nu, (N+1)*sizeof(int));
    c_ptr += (N+1)*sizeof(int);

    (*qp_in)->nb = (int *) c_ptr;
    memcpy(c_ptr, nb, (N+1)*sizeof(int));
    c_ptr += (N+1)*sizeof(int);

    (*qp_in)->nc = (int *) c_ptr;
    memcpy(c_ptr, nc, (N+1)*sizeof(int));
    c_ptr += (N+1)*sizeof(int);

    // assign double pointers
    (*qp_in)->A = (double **) c_ptr;
    c_ptr += N*sizeof(double *);

    (*qp_in)->B = (double **) c_ptr;
    c_ptr += N*sizeof(double *);

    (*qp_in)->b = (double **) c_ptr;
    c_ptr += N*sizeof(double *);

    (*qp_in)->Q = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->S = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->R = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->q = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->r = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->idxb = (int **) c_ptr;
    c_ptr += (N+1)*sizeof(int *);

    (*qp_in)->lb = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->ub = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->Cx = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->Cu = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->lc = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_in)->uc = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    // assign pointers to ints
    for (int k = 0; k < N+1; k++) {
        (*qp_in)->idxb[k] = (int *) c_ptr;
        c_ptr += nb[k]*sizeof(int);
    }

    // align data
    align_char_to(8, &c_ptr);

    // assign pointers to doubles
    c_ptr_QPdata = c_ptr;

    for (int k = 0; k < N+1; k++) {
        // assert((size_t)c_ptr % 8 == 0);

        if (k < N) {
            (*qp_in)->A[k] = (double *) c_ptr;
            c_ptr += nx[k+1]*nx[k]*sizeof(double);

            (*qp_in)->B[k] = (double *) c_ptr;
            c_ptr += nx[k+1]*nu[k]*sizeof(double);

            (*qp_in)->b[k] = (double *) c_ptr;
            c_ptr += nx[k+1]*sizeof(double);
        }

        (*qp_in)->Q[k] = (double *) c_ptr;
        c_ptr += nx[k]*nx[k]*sizeof(double);

        (*qp_in)->S[k] = (double *) c_ptr;
        c_ptr += nu[k]*nx[k]*sizeof(double);

        (*qp_in)->R[k] = (double *) c_ptr;
        c_ptr += nu[k]*nu[k]*sizeof(double);

        (*qp_in)->q[k] = (double *) c_ptr;
        c_ptr += nx[k]*sizeof(double);

        (*qp_in)->r[k] = (double *) c_ptr;
        c_ptr += nu[k]*sizeof(double);

        (*qp_in)->lb[k] = (double *) c_ptr;
        c_ptr += nb[k]*sizeof(double);

        (*qp_in)->ub[k] = (double *) c_ptr;
        c_ptr += nb[k]*sizeof(double);

        (*qp_in)->Cx[k] = (double *) c_ptr;
        c_ptr += nc[k]*nx[k]*sizeof(double);

        (*qp_in)->Cu[k] = (double *) c_ptr;
        c_ptr += nc[k]*nu[k]*sizeof(double);

        (*qp_in)->lc[k] = (double *) c_ptr;
        c_ptr += nc[k]*sizeof(double);

        (*qp_in)->uc[k] = (double *) c_ptr;
        c_ptr += nc[k]*sizeof(double);
    }

    // set QP data to zero (mainly for valgrind)
    for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
    *idx = 0;

    return c_ptr;
}



int col_maj_ocp_qp_out_calculate_size(ocp_qp_dims *dims)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *nc = dims->ng;

    int size = sizeof(col_maj_ocp_qp_out);

    size += 3*(N+1)*sizeof(double *);  // u, x, lam
    size += N*sizeof(double *);  // pi

    for (int k = 0; k < N+1; k++) {
        size += (nx[k] + nu[k])*sizeof(double);  // u, x
        if (k < N)
            size += (nx[k+1])*sizeof(double);  // pi
        size += 2*(nb[k] + nc[k])*sizeof(double);  // lam
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *col_maj_ocp_qp_out_assign(ocp_qp_dims *dims, col_maj_ocp_qp_out **qp_out, void *ptr)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *nc = dims->ng;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_out = (col_maj_ocp_qp_out *) c_ptr;
    c_ptr += sizeof(col_maj_ocp_qp_out);

    // assign double pointers
    (*qp_out)->x = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_out)->u = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    (*qp_out)->pi = (double **) c_ptr;
    c_ptr += N*sizeof(double *);

    (*qp_out)->lam = (double **) c_ptr;
    c_ptr += (N+1)*sizeof(double *);

    // align data
    align_char_to(8, &c_ptr);

    // NOTE(dimitris): splitted the loops below to be able to print primal/dual solution at once

    // assign pointers to QP solution
    for (int k = 0; k < N+1; k++) {

        (*qp_out)->x[k] = (double *) c_ptr;
        c_ptr += nx[k]*sizeof(double);

        (*qp_out)->u[k] = (double *) c_ptr;
        c_ptr += nu[k]*sizeof(double);
    }

    for (int k = 0; k < N; k++) {
        (*qp_out)->pi[k] = (double *) c_ptr;
        c_ptr += nx[k+1]*sizeof(double);
    }

    for (int k = 0; k < N+1; k++) {
        (*qp_out)->lam[k] = (double *) c_ptr;
        c_ptr += 2*(nb[k] + nc[k])*sizeof(double);
    }
    return c_ptr;
}



void convert_from_col_maj_ocp_qp_in(ocp_qp_dims *dims, col_maj_ocp_qp_in *cm_qp_in, ocp_qp_in *qp_in)
{
    qp_in->size->N = cm_qp_in->N;

    // bounds and idxb
    for (int ii = 0; ii <= dims->N; ii++)
    {
        qp_in->size->nbx[ii] = 0;
        qp_in->size->nbu[ii] = 0;
        for (int jj = 0; jj < cm_qp_in->nb[ii]; jj++)
        {

            if (cm_qp_in->idxb[ii][jj] < cm_qp_in->nx[ii])  // state bound
            {
                qp_in->size->nbx[ii]++;
                qp_in->idxb[ii][jj] = cm_qp_in->idxb[ii][jj] + cm_qp_in->nu[ii];
            } else
            {
                qp_in->size->nbu[ii]++;
                qp_in->idxb[ii][jj] = cm_qp_in->idxb[ii][jj] - cm_qp_in->nx[ii];
            }
        }
        qp_in->size->nb[ii] = qp_in->size->nbx[ii] + qp_in->size->nbu[ii];
        assert(qp_in->size->nb[ii] == cm_qp_in->nb[ii]);
    }


    for (int ii = 0; ii <= dims->N; ii++)
    {
        // rest of dimensions
        qp_in->size->nx[ii] = cm_qp_in->nx[ii];
        qp_in->size->nu[ii] = cm_qp_in->nu[ii];
        qp_in->size->ng[ii] = cm_qp_in->nc[ii];
        qp_in->size->ns[ii] = 0;

        // objective

        // dynamics
        if (ii < dims->N)
        {

        }
    }

}



void convert_to_col_maj_ocp_qp_out(ocp_qp_dims *dims, ocp_qp_out *qp_out, col_maj_ocp_qp_out *cm_qp_out)
{
    for (int ii = 0; ii <= dims->N; ii++)
    {
		d_cvt_strvec2vec(dims->nu[ii], &qp_out->ux[ii], 0, cm_qp_out->u[ii]);
        d_cvt_strvec2vec(dims->nx[ii], &qp_out->ux[ii], dims->nu[ii], cm_qp_out->x[ii]);

        if (ii < dims->N)
        {
            d_cvt_strvec2vec(dims->nx[ii+1], &qp_out->pi[ii], 0, cm_qp_out->pi[ii]);
        }

        // TODO(dimitris): change to new convention for the col_maj interface
        d_cvt_strvec2vec(2*dims->nb[ii]+2*dims->ng[ii], &qp_out->lam[ii], 0, cm_qp_out->lam[ii]);
    }

    // col_maj_ocp_qp_out *sol = cm_qp_out;

    // // dummy qp_in
    // ocp_qp_in qp_in;
    // qp_in.N = dims->N;
    // qp_in.nx = dims->nx;
    // qp_in.nu = dims->nu;
    // qp_in.nb = dims->nb;
    // qp_in.ng = dims->ng;
    // qp_in.ns = dims->ns;

    // d_cvt_ocp_qp_sol_to_colmaj(&qp_in, qp_out, sol->u, sol->x, ls, us, sol->pi, lam_lb, lam_ub, lam_lg, lam_ug, lam_ls, lam_us);
}

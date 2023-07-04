/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


// external
// TODO(dimitris): remove memcpy to avoid this dependency?
#include <assert.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"

// hpipm
// #include "hpipm/include/hpipm_d_ocp_qp.h"
// #include "hpipm/include/hpipm_d_ocp_qp_sol.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/utils/mem.h"
// #include "utils/types.h"
// #include "ocp_qp/ocp_qp_hpipm.h"



acados_size_t colmaj_ocp_qp_in_calculate_size(ocp_qp_dims *dims)
{
    size_t N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *nc = dims->ng;

    acados_size_t size = sizeof(colmaj_ocp_qp_in);

    size += 4 * (N + 1) * sizeof(int);        // nx, nu, nb, nc
    size += 3 * N * sizeof(double *);         // A, B, b
    size += 11 * (N + 1) * sizeof(double *);  // ...
    size += 1 * (N + 1) * sizeof(int *);      // idxb

    for (int k = 0; k < N + 1; k++)
    {
        if (k < N)
        {
            size += nx[k + 1] * nx[k] * sizeof(double);  // A
            size += nx[k + 1] * nu[k] * sizeof(double);  // B
            size += nx[k + 1] * sizeof(double);          // b
        }

        size += nx[k] * nx[k] * sizeof(double);  // Q
        size += nu[k] * nx[k] * sizeof(double);  // S
        size += nu[k] * nu[k] * sizeof(double);  // R
        size += nx[k] * sizeof(double);          // q
        size += nu[k] * sizeof(double);          // r
        size += nb[k] * sizeof(int);             // idxb
        size += 2 * nb[k] * sizeof(double);      // lb, ub
        size += nc[k] * nx[k] * sizeof(double);  // Cx
        size += nc[k] * nu[k] * sizeof(double);  // Cu
        size += 2 * nc[k] * sizeof(double);      // lc, uc
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *assign_colmaj_ocp_qp_in(ocp_qp_dims *dims, colmaj_ocp_qp_in **qp_in, void *ptr)
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

    *qp_in = (colmaj_ocp_qp_in *) c_ptr;
    c_ptr += sizeof(colmaj_ocp_qp_in);

    // copy dimensions to workspace
    (*qp_in)->N = N;

    (*qp_in)->nx = (int *) c_ptr;
    memcpy(c_ptr, nx, (N + 1) * sizeof(int));
    c_ptr += (N + 1) * sizeof(int);

    (*qp_in)->nu = (int *) c_ptr;
    memcpy(c_ptr, nu, (N + 1) * sizeof(int));
    c_ptr += (N + 1) * sizeof(int);

    (*qp_in)->nb = (int *) c_ptr;
    memcpy(c_ptr, nb, (N + 1) * sizeof(int));
    c_ptr += (N + 1) * sizeof(int);

    (*qp_in)->nc = (int *) c_ptr;
    memcpy(c_ptr, nc, (N + 1) * sizeof(int));
    c_ptr += (N + 1) * sizeof(int);

    // assign double pointers
    (*qp_in)->A = (double **) c_ptr;
    c_ptr += N * sizeof(double *);

    (*qp_in)->B = (double **) c_ptr;
    c_ptr += N * sizeof(double *);

    (*qp_in)->b = (double **) c_ptr;
    c_ptr += N * sizeof(double *);

    (*qp_in)->Q = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->S = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->R = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->q = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->r = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->idxb = (int **) c_ptr;
    c_ptr += (N + 1) * sizeof(int *);

    (*qp_in)->lb = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->ub = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->Cx = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->Cu = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->lc = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_in)->uc = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    // assign pointers to ints
    for (int k = 0; k < N + 1; k++)
    {
        (*qp_in)->idxb[k] = (int *) c_ptr;
        c_ptr += nb[k] * sizeof(int);
    }

    // align data
    align_char_to(8, &c_ptr);

    // assign pointers to doubles
    c_ptr_QPdata = c_ptr;

    for (int k = 0; k < N + 1; k++)
    {
        // assert((size_t)c_ptr % 8 == 0);

        if (k < N)
        {
            (*qp_in)->A[k] = (double *) c_ptr;
            c_ptr += nx[k + 1] * nx[k] * sizeof(double);

            (*qp_in)->B[k] = (double *) c_ptr;
            c_ptr += nx[k + 1] * nu[k] * sizeof(double);

            (*qp_in)->b[k] = (double *) c_ptr;
            c_ptr += nx[k + 1] * sizeof(double);
        }

        (*qp_in)->Q[k] = (double *) c_ptr;
        c_ptr += nx[k] * nx[k] * sizeof(double);

        (*qp_in)->S[k] = (double *) c_ptr;
        c_ptr += nu[k] * nx[k] * sizeof(double);

        (*qp_in)->R[k] = (double *) c_ptr;
        c_ptr += nu[k] * nu[k] * sizeof(double);

        (*qp_in)->q[k] = (double *) c_ptr;
        c_ptr += nx[k] * sizeof(double);

        (*qp_in)->r[k] = (double *) c_ptr;
        c_ptr += nu[k] * sizeof(double);

        (*qp_in)->lb[k] = (double *) c_ptr;
        c_ptr += nb[k] * sizeof(double);

        (*qp_in)->ub[k] = (double *) c_ptr;
        c_ptr += nb[k] * sizeof(double);

        (*qp_in)->Cx[k] = (double *) c_ptr;
        c_ptr += nc[k] * nx[k] * sizeof(double);

        (*qp_in)->Cu[k] = (double *) c_ptr;
        c_ptr += nc[k] * nu[k] * sizeof(double);

        (*qp_in)->lc[k] = (double *) c_ptr;
        c_ptr += nc[k] * sizeof(double);

        (*qp_in)->uc[k] = (double *) c_ptr;
        c_ptr += nc[k] * sizeof(double);
    }

    // set QP data to zero (mainly for valgrind)
    for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++) *idx = 0;

    return c_ptr;
}



acados_size_t colmaj_ocp_qp_out_calculate_size(ocp_qp_dims *dims)
{
    size_t N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *nc = dims->ng;

    acados_size_t size = sizeof(colmaj_ocp_qp_out);

    size += 3 * (N + 1) * sizeof(double *);  // u, x, lam
    size += N * sizeof(double *);            // pi

    for (size_t k = 0; k < N + 1; k++)
    {
        size += (nx[k] + nu[k]) * sizeof(double);         // u, x
        if (k < N) size += (nx[k + 1]) * sizeof(double);  // pi
        size += 2 * (nb[k] + nc[k]) * sizeof(double);     // lam
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *assign_colmaj_ocp_qp_out(ocp_qp_dims *dims, colmaj_ocp_qp_out **qp_out, void *ptr)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *nc = dims->ng;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_out = (colmaj_ocp_qp_out *) c_ptr;
    c_ptr += sizeof(colmaj_ocp_qp_out);

    // assign double pointers
    (*qp_out)->x = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_out)->u = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    (*qp_out)->pi = (double **) c_ptr;
    c_ptr += N * sizeof(double *);

    (*qp_out)->lam = (double **) c_ptr;
    c_ptr += (N + 1) * sizeof(double *);

    // align data
    align_char_to(8, &c_ptr);

    // NOTE(dimitris): splitted the loops below to be able to print primal/dual solution at once

    // assign pointers to QP solution
    for (int k = 0; k < N + 1; k++)
    {
        (*qp_out)->x[k] = (double *) c_ptr;
        c_ptr += nx[k] * sizeof(double);

        (*qp_out)->u[k] = (double *) c_ptr;
        c_ptr += nu[k] * sizeof(double);
    }

    for (int k = 0; k < N; k++)
    {
        (*qp_out)->pi[k] = (double *) c_ptr;
        c_ptr += nx[k + 1] * sizeof(double);
    }

    for (int k = 0; k < N + 1; k++)
    {
        (*qp_out)->lam[k] = (double *) c_ptr;
        c_ptr += 2 * (nb[k] + nc[k]) * sizeof(double);
    }
    return c_ptr;
}



acados_size_t colmaj_ocp_qp_res_calculate_size(ocp_qp_dims *dims)
{
    size_t N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    acados_size_t size = sizeof(colmaj_ocp_qp_res);

    size += 1 * N * sizeof(double *);         // res_b
    size += 16 * (N + 1) * sizeof(double *);  // everything else

    for (size_t k = 0; k <= N; k++)
    {
        size += nu[k] * sizeof(double);      // res_r
        size += nx[k] * sizeof(double);      // res_q
        size += 2 * ns[k] * sizeof(double);  // res_ls, res_us

        if (k < N) size += nx[k + 1] * sizeof(double);  // res_b

        size += 4 * nb[k] * sizeof(double);  // res_d_lb, res_d_ub, res_m_lb, res_m_ub
        size += 4 * ng[k] * sizeof(double);  // res_d_lg, res_d_ug, res_m_lg, res_m_ug
        size += 4 * ns[k] * sizeof(double);  // res_d_ls, res_d_us, res_m_ls, res_m_us
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



char *assign_colmaj_ocp_qp_res(ocp_qp_dims *dims, colmaj_ocp_qp_res **qp_res, void *ptr)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_res = (colmaj_ocp_qp_res *) c_ptr;
    c_ptr += sizeof(colmaj_ocp_qp_res);

    // assign double pointers
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_r, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_q, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_ls, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_us, &c_ptr);
    assign_and_advance_double_ptrs(N, &(*qp_res)->res_b, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_d_lb, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_d_ub, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_d_lg, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_d_ug, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_d_ls, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_d_us, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_m_lb, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_m_ub, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_m_lg, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_m_ug, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_m_ls, &c_ptr);
    assign_and_advance_double_ptrs(N + 1, &(*qp_res)->res_m_us, &c_ptr);

    // align data
    align_char_to(8, &c_ptr);

    // assign pointers to QP solution
    for (int k = 0; k < N + 1; k++)
    {
        assign_and_advance_double(nu[k], &(*qp_res)->res_r[k], &c_ptr);
        assign_and_advance_double(nx[k], &(*qp_res)->res_q[k], &c_ptr);
        assign_and_advance_double(ns[k], &(*qp_res)->res_ls[k], &c_ptr);
        assign_and_advance_double(ns[k], &(*qp_res)->res_us[k], &c_ptr);

        if (k < N) assign_and_advance_double(nx[k + 1], &(*qp_res)->res_b[k], &c_ptr);

        assign_and_advance_double(nb[k], &(*qp_res)->res_d_lb[k], &c_ptr);
        assign_and_advance_double(nb[k], &(*qp_res)->res_d_ub[k], &c_ptr);
        assign_and_advance_double(ng[k], &(*qp_res)->res_d_lg[k], &c_ptr);
        assign_and_advance_double(ng[k], &(*qp_res)->res_d_ug[k], &c_ptr);
        assign_and_advance_double(ns[k], &(*qp_res)->res_d_ls[k], &c_ptr);
        assign_and_advance_double(ns[k], &(*qp_res)->res_d_us[k], &c_ptr);

        assign_and_advance_double(nb[k], &(*qp_res)->res_m_lb[k], &c_ptr);
        assign_and_advance_double(nb[k], &(*qp_res)->res_m_ub[k], &c_ptr);
        assign_and_advance_double(ng[k], &(*qp_res)->res_m_lg[k], &c_ptr);
        assign_and_advance_double(ng[k], &(*qp_res)->res_m_ug[k], &c_ptr);
        assign_and_advance_double(ns[k], &(*qp_res)->res_m_ls[k], &c_ptr);
        assign_and_advance_double(ns[k], &(*qp_res)->res_m_us[k], &c_ptr);
    }

    return c_ptr;
}



void convert_colmaj_to_ocp_qp_in(colmaj_ocp_qp_in *cm_qp_in, ocp_qp_in *qp_in)
{
    int N = cm_qp_in->N;

    // bounds and idxb
    for (int ii = 0; ii < N + 1; ii++)
    {
        qp_in->dim->nbx[ii] = 0;
        qp_in->dim->nbu[ii] = 0;
        for (int jj = 0; jj < cm_qp_in->nb[ii]; jj++)
        {
            if (cm_qp_in->idxb[ii][jj] < cm_qp_in->nx[ii])
            {  // state bound
                qp_in->dim->nbx[ii]++;
                qp_in->idxb[ii][jj] = cm_qp_in->idxb[ii][jj] + cm_qp_in->nu[ii];
            }
            else
            {
                qp_in->dim->nbu[ii]++;
                qp_in->idxb[ii][jj] = cm_qp_in->idxb[ii][jj] - cm_qp_in->nx[ii];
            }
        }
        qp_in->dim->nb[ii] = qp_in->dim->nbx[ii] + qp_in->dim->nbu[ii];
        assert(qp_in->dim->nb[ii] == cm_qp_in->nb[ii]);
    }

    // rest of dimensions
    qp_in->dim->N = N;

    for (int ii = 0; ii < N + 1; ii++)
    {
        qp_in->dim->nx[ii] = cm_qp_in->nx[ii];
        qp_in->dim->nu[ii] = cm_qp_in->nu[ii];
        qp_in->dim->ng[ii] = cm_qp_in->nc[ii];
        qp_in->dim->ns[ii] = 0;
    }

    // convert to backend qp
    int ii;
    for(ii=0; ii<N; ii++)
    {
        d_ocp_qp_set_A(ii, cm_qp_in->A[ii], qp_in);
        d_ocp_qp_set_B(ii, cm_qp_in->B[ii], qp_in);
        d_ocp_qp_set_b(ii, cm_qp_in->b[ii], qp_in);
        d_ocp_qp_set_R(ii, cm_qp_in->R[ii], qp_in);
        d_ocp_qp_set_S(ii, cm_qp_in->S[ii], qp_in);
        d_ocp_qp_set_Q(ii, cm_qp_in->Q[ii], qp_in);
        d_ocp_qp_set_r(ii, cm_qp_in->r[ii], qp_in);
        d_ocp_qp_set_q(ii, cm_qp_in->q[ii], qp_in);
        d_ocp_qp_set_idxb(ii, cm_qp_in->idxb[ii], qp_in);
        d_ocp_qp_set_lb(ii, cm_qp_in->lb[ii], qp_in);
        d_ocp_qp_set_ub(ii, cm_qp_in->ub[ii], qp_in);
        d_ocp_qp_set_C(ii, cm_qp_in->Cx[ii], qp_in);
        d_ocp_qp_set_D(ii, cm_qp_in->Cu[ii], qp_in);
        d_ocp_qp_set_lg(ii, cm_qp_in->lc[ii], qp_in);
        d_ocp_qp_set_ug(ii, cm_qp_in->uc[ii], qp_in);
    }
    d_ocp_qp_set_R(ii, cm_qp_in->R[ii], qp_in);
    d_ocp_qp_set_S(ii, cm_qp_in->S[ii], qp_in);
    d_ocp_qp_set_Q(ii, cm_qp_in->Q[ii], qp_in);
    d_ocp_qp_set_r(ii, cm_qp_in->r[ii], qp_in);
    d_ocp_qp_set_q(ii, cm_qp_in->q[ii], qp_in);
    d_ocp_qp_set_idxb(ii, cm_qp_in->idxb[ii], qp_in);
    d_ocp_qp_set_lb(ii, cm_qp_in->lb[ii], qp_in);
    d_ocp_qp_set_ub(ii, cm_qp_in->ub[ii], qp_in);
    d_ocp_qp_set_C(ii, cm_qp_in->Cx[ii], qp_in);
    d_ocp_qp_set_D(ii, cm_qp_in->Cu[ii], qp_in);
    d_ocp_qp_set_lg(ii, cm_qp_in->lc[ii], qp_in);
    d_ocp_qp_set_ug(ii, cm_qp_in->uc[ii], qp_in);

}



void convert_ocp_qp_out_to_colmaj(ocp_qp_out *qp_out, colmaj_ocp_qp_out *cm_qp_out)
{
    ocp_qp_dims *dims = qp_out->dim;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        blasfeo_unpack_dvec(dims->nu[ii], &qp_out->ux[ii], 0, cm_qp_out->u[ii], 1);
        blasfeo_unpack_dvec(dims->nx[ii], &qp_out->ux[ii], dims->nu[ii], cm_qp_out->x[ii], 1);

        if (ii < dims->N)
            blasfeo_unpack_dvec(dims->nx[ii + 1], &qp_out->pi[ii], 0, cm_qp_out->pi[ii], 1);

        // TODO(dimitris): change to new convention for the colmaj interface
        blasfeo_unpack_dvec(2 * dims->nb[ii] + 2 * dims->ng[ii], &qp_out->lam[ii], 0,
                            cm_qp_out->lam[ii], 1);
    }

    // colmaj_ocp_qp_out *sol = cm_qp_out;

    // // dummy qp_in
    // ocp_qp_in qp_in;
    // qp_in.N = dims->N;
    // qp_in.nx = dims->nx;
    // qp_in.nu = dims->nu;
    // qp_in.nb = dims->nb;
    // qp_in.ng = dims->ng;
    // qp_in.ns = dims->ns;

    // d_cvt_ocp_qp_sol_to_colmaj(&qp_in, qp_out, sol->u, sol->x, ls, us, sol->pi, lam_lb, lam_ub,
    //                            lam_lg, lam_ug, lam_ls, lam_us);
}



void convert_ocp_qp_res_to_colmaj(ocp_qp_res *qp_res, colmaj_ocp_qp_res *cm_qp_res)
{
    double **res_r = cm_qp_res->res_r;
    double **res_q = cm_qp_res->res_q;
    double **res_ls = cm_qp_res->res_ls;
    double **res_us = cm_qp_res->res_us;
    double **res_b = cm_qp_res->res_b;
    double **res_d_lb = cm_qp_res->res_d_lb;
    double **res_d_ub = cm_qp_res->res_d_ub;
    double **res_d_lg = cm_qp_res->res_d_lg;
    double **res_d_ug = cm_qp_res->res_d_ug;
    double **res_d_ls = cm_qp_res->res_d_ls;
    double **res_d_us = cm_qp_res->res_d_us;
    double **res_m_lb = cm_qp_res->res_m_lb;
    double **res_m_ub = cm_qp_res->res_m_ub;
    double **res_m_lg = cm_qp_res->res_m_lg;
    double **res_m_ug = cm_qp_res->res_m_ug;
    double **res_m_ls = cm_qp_res->res_m_ls;
    double **res_m_us = cm_qp_res->res_m_us;

    ocp_qp_res_compute_nrm_inf(qp_res, cm_qp_res->res_nrm_inf);

    d_ocp_qp_res_get_all(qp_res, res_r, res_q, res_ls, res_us, res_b, res_d_lb, res_d_ub,
                               res_d_lg, res_d_ug, res_d_ls, res_d_us, res_m_lb, res_m_ub, res_m_lg,
                               res_m_ug, res_m_ls, res_m_us);
}

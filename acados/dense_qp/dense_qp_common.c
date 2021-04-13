/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stddef.h>
#include <string.h>

// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_kkt.h"
#include "hpipm/include/hpipm_d_dense_qp_res.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp_dim.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"




/************************************************
 * config
 ************************************************/

acados_size_t dense_qp_solver_config_calculate_size()
{
    acados_size_t size = 0;

    size += sizeof(qp_solver_config);

    return size;
}



qp_solver_config *dense_qp_solver_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    qp_solver_config *config = (qp_solver_config *) c_ptr;
    c_ptr += sizeof(qp_solver_config);

    return config;
}



/************************************************
 * dims
 ************************************************/

acados_size_t dense_qp_dims_calculate_size()
{
    acados_size_t size = sizeof(dense_qp_dims);

    size += d_dense_qp_dim_memsize();

    return size;
}



dense_qp_dims *dense_qp_dims_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_dims *dims = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    d_dense_qp_dim_create(dims, c_ptr);
    c_ptr += d_dense_qp_dim_memsize();

    assert((char *) raw_memory + dense_qp_dims_calculate_size() == c_ptr);

    return dims;
}


void dense_qp_dims_set(void *config_, void *dims_, const char *field, const int* value)
{
    // wrap hpipm function
    // char field_copy[MAX_STR_LEN];
    // strcpy(field_copy, field);
    char *field_copy = (char *) field;

    dense_qp_dims *dims = (dense_qp_dims *) dims_;

    d_dense_qp_dim_set(field_copy, *value, dims);
}



/************************************************
 * in
 ************************************************/

acados_size_t dense_qp_in_calculate_size(dense_qp_dims *dims)
{
    acados_size_t size = sizeof(dense_qp_in);
    size += sizeof(dense_qp_dims);
    size += d_dense_qp_memsize(dims);
    size += 8;  // align for d_dense_qp
    make_int_multiple_of(8, &size);

    return size;
}



dense_qp_in *dense_qp_in_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_in *qp_in = (dense_qp_in *) c_ptr;
    c_ptr += sizeof(dense_qp_in);

    align_char_to(8, &c_ptr);
    d_dense_qp_create(dims, qp_in, c_ptr);
    c_ptr += d_dense_qp_memsize(dims);

    qp_in->dim = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    qp_in->dim->nv = dims->nv;
    qp_in->dim->ne = dims->ne;
    qp_in->dim->nb = dims->nb;
    qp_in->dim->ng = dims->ng;
    qp_in->dim->ns = dims->ns;
    qp_in->dim->nsb = dims->nsb;
    qp_in->dim->nsg = dims->nsg;

    assert((char *) raw_memory + dense_qp_in_calculate_size(dims) >= c_ptr);

    return qp_in;
}



/************************************************
 * out
 ************************************************/

acados_size_t dense_qp_out_calculate_size(dense_qp_dims *dims)
{
    acados_size_t size = sizeof(dense_qp_out);
    size += d_dense_qp_sol_memsize(dims);
    size += sizeof(qp_info);
    size += 8;
    make_int_multiple_of(8, &size);

    return size;
}



dense_qp_out *dense_qp_out_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_out *qp_out = (dense_qp_out *) c_ptr;
    c_ptr += sizeof(dense_qp_out);
    align_char_to(8, &c_ptr);

    d_dense_qp_sol_create(dims, qp_out, c_ptr);
    c_ptr += d_dense_qp_sol_memsize(dims);

    qp_out->misc = (void *) c_ptr;
    c_ptr += sizeof(qp_info);

    assert((char *) raw_memory + dense_qp_out_calculate_size(dims) >= c_ptr);

    return qp_out;
}



void dense_qp_out_get(dense_qp_out *out, const char *field, void *value)
{
    if(!strcmp(field, "qp_info"))
    {
        qp_info **ptr = value;
        *ptr = out->misc;
    }
    else
    {
        printf("\nerror: dense_qp_out_get: field %s not available\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * res
 ************************************************/

acados_size_t dense_qp_res_calculate_size(dense_qp_dims *dims)
{
    acados_size_t size = sizeof(dense_qp_res);
    size += d_dense_qp_res_memsize(dims);

    make_int_multiple_of(8, &size);
    return size;
}



dense_qp_res *dense_qp_res_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res *qp_res = (dense_qp_res *) c_ptr;
    c_ptr += sizeof(dense_qp_res);

    align_char_to(8, &c_ptr);
    d_dense_qp_res_create(dims, qp_res, c_ptr);
    c_ptr += d_dense_qp_res_memsize(dims);
    assert((char *) raw_memory + dense_qp_res_calculate_size(dims) >= c_ptr);

    return qp_res;
}



acados_size_t dense_qp_res_workspace_calculate_size(dense_qp_dims *dims)
{
    acados_size_t size = sizeof(dense_qp_res_ws);

    size += d_dense_qp_res_ws_memsize(dims);
    make_int_multiple_of(8, &size);

    return size;
}



dense_qp_res_ws *dense_qp_res_workspace_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res_ws *res_ws = (dense_qp_res_ws *) c_ptr;
    c_ptr += sizeof(dense_qp_res_ws);

    align_char_to(8, &c_ptr);

    d_dense_qp_res_ws_create(dims, res_ws, c_ptr);
    c_ptr += d_dense_qp_res_ws_memsize(dims);

    assert((char *) raw_memory + dense_qp_res_workspace_calculate_size(dims) >= c_ptr);

    return res_ws;
}



void dense_qp_compute_t(dense_qp_in *qp_in, dense_qp_out *qp_out)
{
    int nvd = qp_in->dim->nv;
    // int ned = qp_in->dim->ne;
    int nbd = qp_in->dim->nb;
    int ngd = qp_in->dim->ng;
    int nsd = qp_in->dim->ns;

    int *idxb = qp_in->idxb;
    int *idxs_rev = qp_in->idxs_rev;

    int ii, idx;

    // compute slacks for bounds
    blasfeo_dvecex_sp(nbd, 1.0, idxb, qp_out->v, 0, qp_out->t, nbd+ngd);
    blasfeo_daxpby(nbd, 1.0, qp_out->t, nbd+ngd, -1.0, qp_in->d, 0, qp_out->t, 0);
    blasfeo_daxpby(nbd, -1.0, qp_out->t, nbd+ngd, -1.0, qp_in->d, nbd + ngd, qp_out->t, nbd + ngd);

    // compute slacks for general constraints
    blasfeo_dgemv_t(nvd, ngd, 1.0, qp_in->Ct, 0, 0, qp_out->v, 0, -1.0, qp_in->d, nbd, qp_out->t,
                    nbd);
    blasfeo_dgemv_t(nvd, ngd, -1.0, qp_in->Ct, 0, 0, qp_out->v, 0, -1.0, qp_in->d, 2 * nbd + ngd,
                    qp_out->t, 2 * nbd + ngd);

    // soft
//    blasfeo_dvecad_sp(nsd, 1.0, qp_out->v, nvd, idxs, qp_out->t, 0);
//    blasfeo_dvecad_sp(nsd, 1.0, qp_out->v, nvd+nsd, idxs, qp_out->t, nbd+ngd);
    for(ii=0; ii<nbd+ngd; ii++)
    {
        idx = idxs_rev[ii];
        if(idx!=-1)
        {
        BLASFEO_DVECEL(qp_out->t, ii)         += BLASFEO_DVECEL(qp_out->v, nvd+idx);
        BLASFEO_DVECEL(qp_out->t, nbd+ngd+ii) += BLASFEO_DVECEL(qp_out->v, nvd+nsd+idx);
        }
    }
    blasfeo_dveccp(2*nsd, qp_out->v, nvd, qp_out->t, 2*nbd+2*ngd);
    blasfeo_daxpy(2*nsd, -1.0, qp_in->d, 2*nbd+2*ngd, qp_out->t, 2*nbd+2*ngd, qp_out->t,
        2*nbd+2*ngd);

}



void dense_qp_res_compute(dense_qp_in *qp_in, dense_qp_out *qp_out, dense_qp_res *qp_res,
                          dense_qp_res_ws *res_ws)
{
    qp_info *info = (qp_info *) qp_out->misc;

    if (info->t_computed == 0)
    {
        dense_qp_compute_t(qp_in, qp_out);
        info->t_computed = 1;
    }

    // compute residuals
    d_dense_qp_res_compute(qp_in, qp_out, qp_res, res_ws);
}



void dense_qp_res_compute_nrm_inf(dense_qp_res *qp_res, double res[4])
{
    int nv = qp_res->dim->nv;
    int nb = qp_res->dim->nb;
    int ne = qp_res->dim->ne;
    int ng = qp_res->dim->ng;
    int ns = qp_res->dim->ns;

    blasfeo_dvecnrm_inf(nv + 2 * ns, qp_res->res_g, 0, &res[0]);
    blasfeo_dvecnrm_inf(ne, qp_res->res_b, 0, &res[1]);
    blasfeo_dvecnrm_inf(2 * nb + 2 * ng + 2 * ns, qp_res->res_d, 0, &res[2]);
    blasfeo_dvecnrm_inf(2 * nb + 2 * ng + 2 * ns, qp_res->res_m, 0, &res[3]);
}



void dense_qp_stack_slacks_dims(dense_qp_dims *in, dense_qp_dims *out)
{
    out->nv = in->nv + 2 * in->ns;
    out->ne = in->ne;
    out->nb = in->nb - in->nsb + 2 * in->ns;
    out->ng = in->ns > 0 ? in->ng + in->nsb : in->ng;
    out->ns = 0;
    out->nsb = 0;
    out->nsg = 0;
}



void dense_qp_stack_slacks(dense_qp_in *in, dense_qp_in *out)
{
    int nv = in->dim->nv;
    int ne = in->dim->ne;
    int nb = in->dim->nb;
    int ng = in->dim->ng;
    int ns = in->dim->ns;
    int nsb = in->dim->nsb;
    // int nsg = in->dim->nsg;
    int *idxs_rev = in->idxs_rev;
    int *idxb = in->idxb;

    int nv2 = out->dim->nv;
    int ne2 = out->dim->ne;
    int nb2 = out->dim->nb;
    int ng2 = out->dim->ng;

    assert(nv2 == nv+2*ns && "Dimensions are wrong!");
    assert(nb2 == nb-nsb+2*ns && "Dimensions are wrong!");
    assert(ng2 == ng+nsb && "Dimensions are wrong!");

    // set matrices to 0.0
    blasfeo_dgese(nv2, nv2, 0.0, out->Hv, 0, 0);
    blasfeo_dgese(ne2, nv2, 0.0, out->A, 0, 0);
    blasfeo_dgese(nv2, ng2, 0.0, out->Ct, 0, 0);

    // copy in->Hv to upper left corner of out->Hv, out->Hv = [in->Hv 0; 0 0]
    blasfeo_dgecp(nv, nv, in->Hv, 0, 0, out->Hv, 0, 0);

    // copy in->Z to main diagonal of out->Hv, out->Hv = [in->Hv 0; 0 Z]
    blasfeo_ddiain(2 * ns, 1.0, in->Z, 0, out->Hv, nv, nv);

    // copy in->gz to out->gz
    blasfeo_dveccp(nv + 2 * ns, in->gz, 0, out->gz, 0);

    if (ne > 0)
    {
        // copy in->A to out->A
        blasfeo_dgecp(ne, nv, in->A, 0, 0, out->A, 0, 0);

        // copy in->b to out->b
        blasfeo_dveccp(ne, in->b, 0, out->b, 0);
    }

    // copy in->Ct to out->Ct, out->Ct = [in->Ct 0; 0 0]
    blasfeo_dgecp(nv, ng, in->Ct, 0, 0, out->Ct, 0, 0);

    if (ns > 0)
    {
        // set flags for non-softened box constraints
        // use out->m temporarily for this
        for (int ii = 0; ii < nb; ii++)
        {
            // TODO(bnovoselnik): pick up some workspace for this
            BLASFEO_DVECEL(out->m, ii) = 1.0;
        }

        int col_b = ng;
//        for (int ii = 0; ii < ns; ii++)
//        {
//            int js = idxs[ii];
        for(int js=0; js<nb+ng; js++)
        {
            int ii = idxs_rev[js];
            if(ii!=-1)
            {

                // int idx_v_ls0 = nv+ii;
                // int idx_v_us0 = nv+ns+ii;
                int idx_v_ls1 = nv+ii;
                int idx_v_us1 = nv+ns+ii;

                int idx_d_ls0 = js;
                int idx_d_us0 = nb+ng+js;
                int idx_d_ls1;
                int idx_d_us1;

                if (js < nb)  // soft box constraint
                {
                    // index of a soft box constraint
                    int jv = idxb[js];

                    idx_d_ls1 = nb2+col_b;
                    idx_d_us1 = 2*nb2+ng2+col_b;

                    // softened box constraint, set its flag to -1
                    BLASFEO_DVECEL(out->m, js) = -1.0;

                    // insert softened box constraint into out->Ct, lb_i <= x_i + sl_i - su_i <= ub_i
                    BLASFEO_DMATEL(out->Ct, jv, col_b) = 1.0;
                    BLASFEO_DMATEL(out->Ct, idx_v_ls1, col_b) = +1.0;
                    BLASFEO_DMATEL(out->Ct, idx_v_us1, col_b) = -1.0;
                    BLASFEO_DVECEL(out->d, idx_d_ls1) = BLASFEO_DVECEL(in->d, idx_d_ls0);
                    BLASFEO_DVECEL(out->d, idx_d_us1) = -BLASFEO_DVECEL(in->d, idx_d_us0);

                    col_b++;
                }
                else  // soft general constraint
                {
                    // index of a soft general constraint
                    int col_g = js - nb;

                    // soft general constraint, lg_i <= C_i x + sl_i - su_i <= ug_i
                    BLASFEO_DMATEL(out->Ct, idx_v_ls1, col_g) = +1.0;
                    BLASFEO_DMATEL(out->Ct, idx_v_us1, col_g) = -1.0;
                }

                // slack variables have box constraints
                out->idxb[nb-nsb+ii] = ii + nv;
                out->idxb[nb-nsb+ns+ii] = ii + nv + ns;
            }
        }

        int k_nsb = 0;
        for (int ii = 0; ii < nb; ii++)
        {
            if (BLASFEO_DVECEL(out->m, ii) > 0)
            {
                // copy nonsoftened box constraint bounds to out->d
                BLASFEO_DVECEL(out->d, k_nsb) = BLASFEO_DVECEL(in->d, ii);
                BLASFEO_DVECEL(out->d, nb2+ng2+k_nsb) = -BLASFEO_DVECEL(in->d, nb+ng+ii);
                out->idxb[k_nsb] = ii;
                k_nsb++;
            }
        }

        assert(k_nsb == nb-nsb && "Dimensions are wrong!");

        // copy ls and us to out->lb, jump over nonsoftened box constraints
        blasfeo_dveccp(2*ns, in->d, 2*nb+2*ng, out->d, k_nsb);

        // for slack variables out->ub is +INFTY
        blasfeo_dvecse(2*ns, 1.0e6, out->d, nb2+ng2+k_nsb);

        // copy in->lg to out->lg
        blasfeo_dveccp(ng, in->d, nb, out->d, nb2);

        // copy in->ug to out->ug
        blasfeo_dveccpsc(ng, -1.0, in->d, 2*nb+ng, out->d, 2*nb2+ng2);

        // flip signs for ub and ug
        blasfeo_dvecsc(nb2+ng2, -1.0, out->d, nb2+ng2);

        // set out->m to 0.0
        blasfeo_dvecse(2*nb2+2*ng2, 0.0, out->m, 0);
    }
    else
    {
        blasfeo_dveccp(2*nb+2*ng, in->d, 0, out->d, 0);
        blasfeo_dveccp(2*nb+2*ng, in->m, 0, out->m, 0);
        for (int ii = 0; ii < nb; ii++) out->idxb[ii] = in->idxb[ii];
    }
}



void dense_qp_unstack_slacks(dense_qp_out *in, dense_qp_in *qp_out, dense_qp_out *out)
{
    int nv = qp_out->dim->nv;
    int ne = qp_out->dim->ne;
    int nb = qp_out->dim->nb;
    int ng = qp_out->dim->ng;
    int ns = qp_out->dim->ns;
    int nsb = qp_out->dim->nsb;
    // int nsg = qp_out->dim->nsg;

    int *idxs_rev = qp_out->idxs_rev;

    int nv2 = in->dim->nv;
    // int ne2 = in->dim->ne;
    int nb2 = in->dim->nb;
    int ng2 = in->dim->ng;

    UNUSED(nsb);
    UNUSED(nv2);

    assert(nv2 == nv+2*ns && "Dimensions are wrong!");
    assert(nb2 == nb-nsb+2*ns && "Dimensions are wrong!");
    assert(ng2 == 2*ng+2*nsb && "Dimensions are wrong!");

    // inequality constraints multipliers and slacks
    if (ns > 0)
    {

        // set flags for non-softened box constraints
        // use out->v temporarily for this
        // XXX assume that nb<=nv
        assert(nv+2*ns >= nb);

        for (int ii = 0; ii < nb; ii++)
        {
            // TODO(bnovoselnik): pick up some workspace for this
            BLASFEO_DVECEL(out->v, ii) = 1.0;
        }

        int col_b = 2*ng;
//        for (int ii = 0; ii < ns; ii++)
//        {
//            int js = idxs[ii];
        for(int js=0; js<nb+ng; js++)
        {
            int ii = idxs_rev[js];
            if(ii!=-1)
            {

                int idx_d_ls0 = js;
                int idx_d_us0 = nb+ng+js;
                int idx_d_ls1;
                int idx_d_us1;

                if (js < nb)
                {
                    // softened box constraint, set its flag to -1
                    BLASFEO_DVECEL(out->v, js) = -1.0;

                    idx_d_ls1 = nb2+col_b;
                    idx_d_us1 = 2*nb2+ng2+col_b;

                    BLASFEO_DVECEL(out->lam, idx_d_ls0) = BLASFEO_DVECEL(in->lam, idx_d_ls1);
                    BLASFEO_DVECEL(out->lam, idx_d_us0) = BLASFEO_DVECEL(in->lam, idx_d_us1);

                    BLASFEO_DVECEL(out->t, idx_d_ls0) = BLASFEO_DVECEL(in->t, idx_d_ls1);
                    BLASFEO_DVECEL(out->t, idx_d_us0) = BLASFEO_DVECEL(in->t, idx_d_us1);

                    col_b++;
                }
            }
        }

        int k_nsb = 0;  // number of non-softed bounds
        for (int ii = 0; ii < nb; ii++)
        {
            int idx_d_ls0 = ii;
            int idx_d_us0 = nb+ng+ii;
            int idx_d_ls1;
            int idx_d_us1;

            if (BLASFEO_DVECEL(out->v, ii) > 0)
            {
                idx_d_ls1 = k_nsb;
                idx_d_us1 = nb2+ng2+k_nsb;

                BLASFEO_DVECEL(out->lam, idx_d_ls0) = BLASFEO_DVECEL(in->lam, idx_d_ls1);
                BLASFEO_DVECEL(out->lam, idx_d_us0) = BLASFEO_DVECEL(in->lam, idx_d_us1);

                BLASFEO_DVECEL(out->t, idx_d_ls0) = BLASFEO_DVECEL(in->t, idx_d_ls1);
                BLASFEO_DVECEL(out->t, idx_d_us0) = BLASFEO_DVECEL(in->t, idx_d_us1);

                k_nsb++;
            }
        }

        assert(k_nsb == nb-nsb && "Dimensions are wrong!");

        blasfeo_dveccp(2*ns, in->t, k_nsb, out->t, 2*nb+2*ng);
        blasfeo_dveccp(ng, in->t, 2*nb2+ng2, out->t, 2*nb+ng);
        blasfeo_dveccp(ng, in->t, nb2, out->t, nb);

        blasfeo_dveccp(2*ns, in->lam, k_nsb, out->lam, 2*nb+2*ng);
        blasfeo_dveccp(ng, in->lam, 2*nb2+ng2, out->lam, 2*nb+ng);
        blasfeo_dveccp(ng, in->lam, nb2, out->lam, nb);

        // variables
        blasfeo_dveccp(nv+2*ns, in->v, 0, out->v, 0);

        // equality constraints multipliers
        blasfeo_dveccp(ne, in->pi, 0, out->pi, 0);

    }

    return;
}

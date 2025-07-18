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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp_kkt.h"
#include "hpipm/include/hpipm_d_ocp_qp_res.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"


/************************************************
 * config
 ************************************************/

acados_size_t ocp_qp_solver_config_calculate_size()
{
    acados_size_t size = 0;

    size += sizeof(qp_solver_config);

    return size;
}



qp_solver_config *ocp_qp_solver_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    qp_solver_config *config = (qp_solver_config *) c_ptr;
    c_ptr += sizeof(qp_solver_config);

    return config;
}



acados_size_t ocp_qp_condensing_config_calculate_size()
{
    acados_size_t size = 0;

    size += sizeof(ocp_qp_xcond_config);

    return size;
}



ocp_qp_xcond_config *ocp_qp_condensing_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_qp_xcond_config *config = (ocp_qp_xcond_config *) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_config);

    return config;
}



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_qp_dims_calculate_size(int N)
{
    acados_size_t size = sizeof(ocp_qp_dims);

    size += d_ocp_qp_dim_memsize(N);
    make_int_multiple_of(8, &size);

    size += 8;

    return size;
}



ocp_qp_dims *ocp_qp_dims_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_dims *dims = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);
    align_char_to(8, &c_ptr);

    d_ocp_qp_dim_create(N, dims, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);

    for (int i=0; i<=N; i++)
    {
        dims->nsbu[i] = 0;
        dims->nsbx[i] = 0;
        dims->nsg[i] = 0;
        dims->ns[i] = 0;
    }

    assert((char *) raw_memory + ocp_qp_dims_calculate_size(N) >= c_ptr);

    return dims;
}



void ocp_qp_dims_set(void *config_, void *dims, int stage, const char *field, int* value)
{
    char *field_copy = (char *) field;

    d_ocp_qp_dim_set(field_copy, stage, *value, dims);

    return;
}


void ocp_qp_dims_get(void *config_, void *dims, int stage, const char *field, int* value)
{
    char *field_copy = (char *) field;

    d_ocp_qp_dim_get(dims, field_copy, stage, value);

    return;
}


static int ocp_qp_dims_get_ni(ocp_qp_dims *dims, int stage)
{
    return dims->nbu[stage] + dims->nbx[stage] + dims->ng[stage] + dims->ns[stage];
}

/************************************************
 * in
 ************************************************/

acados_size_t ocp_qp_in_calculate_size(ocp_qp_dims *dims)
{
    acados_size_t size = sizeof(ocp_qp_in);
    size += d_ocp_qp_memsize(dims);
    size += ocp_qp_dims_calculate_size(dims->N);  // TODO(all): remove !!!
    size += 2*8; // aligns

    make_int_multiple_of(8, &size);

    return size;
}



ocp_qp_in *ocp_qp_in_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    align_char_to(8, &c_ptr);
    ocp_qp_in *qp_in = (ocp_qp_in *) c_ptr;
    c_ptr += sizeof(ocp_qp_in);
    align_char_to(8, &c_ptr);

    d_ocp_qp_create(dims, qp_in, c_ptr);
    c_ptr += d_ocp_qp_memsize(dims);

    ocp_qp_dims *dims_copy = ocp_qp_dims_assign(dims->N, c_ptr);  // TODO(all): remove !!!
    c_ptr += ocp_qp_dims_calculate_size(dims->N);                 // TODO(all): remove !!!

    d_ocp_qp_dim_copy_all(dims, dims_copy);

    qp_in->dim = dims_copy;

    assert((char *) raw_memory + ocp_qp_in_calculate_size(dims) >= c_ptr);

    return qp_in;
}



/************************************************
 * out
 ************************************************/

acados_size_t ocp_qp_out_calculate_size(ocp_qp_dims *dims)
{
    acados_size_t size = sizeof(ocp_qp_out);
    size += d_ocp_qp_sol_memsize(dims);
    size += ocp_qp_dims_calculate_size(dims->N);  // TODO(all): remove !!!
    size += sizeof(qp_info); // TODO move to memory !!!

    size += 2*8; // align
    make_int_multiple_of(8, &size);

    return size;
}



ocp_qp_out *ocp_qp_out_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    align_char_to(8, &c_ptr);
    ocp_qp_out *qp_out = (ocp_qp_out *) c_ptr;
    c_ptr += sizeof(ocp_qp_out);
    align_char_to(8, &c_ptr);

    d_ocp_qp_sol_create(dims, qp_out, c_ptr);
    c_ptr += d_ocp_qp_sol_memsize(dims);

    qp_out->misc = (void *) c_ptr; // TODO move to memory !!!
    c_ptr += sizeof(qp_info); // TODO move to memory !!!
    align_char_to(8, &c_ptr);

    ocp_qp_dims *dims_copy = ocp_qp_dims_assign(dims->N, c_ptr);  // TODO(all): remove !!!
    c_ptr += ocp_qp_dims_calculate_size(dims->N);                 // TODO(all): remove !!!

    d_ocp_qp_dim_copy_all(dims, dims_copy);

    qp_out->dim = dims_copy;

    assert((char *) raw_memory + ocp_qp_out_calculate_size(dims) >= c_ptr);

    return qp_out;
}


double ocp_qp_out_compute_primal_nrm_inf(ocp_qp_out* qp_out)
{
    double res = 0;
    double res_stage = 0;
    ocp_qp_dims *dims = qp_out->dim;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    for (int i = 0; i <= N; i++)
    {
        blasfeo_dvecnrm_inf(nx[i]+nu[i]+2*ns[i], qp_out->ux+i, 0, &res_stage);
        res = res > res_stage ? res : res_stage;
    }
    return res;
}

/************************************************
 * seed
 ************************************************/

acados_size_t ocp_qp_seed_calculate_size(ocp_qp_dims *dims)
{
    acados_size_t size = sizeof(ocp_qp_seed);
    size += d_ocp_qp_seed_memsize(dims);
    return size;
}



ocp_qp_seed *ocp_qp_seed_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_seed *qp_seed = (ocp_qp_seed *) c_ptr;
    c_ptr += sizeof(ocp_qp_seed);

    d_ocp_qp_seed_create(dims, qp_seed, c_ptr);
    c_ptr += d_ocp_qp_seed_memsize(dims);

    assert((char *) raw_memory + ocp_qp_seed_calculate_size(dims) == c_ptr);

    return qp_seed;
}

double ocp_qp_out_compute_dual_nrm_inf(ocp_qp_out* qp_out)
{
    double res = 0;
    double res_stage = 0;
    ocp_qp_dims *dims = qp_out->dim;
    int N = dims->N;
    int *nx = dims->nx;

    for (int i = 0; i < N; i++)
    {
        blasfeo_dvecnrm_inf(nx[i+1], qp_out->pi+i, 0, &res_stage);
        res = res > res_stage ? res : res_stage;
    }
    for (int i = 0; i <= N; i++)
    {
        int ni_stage = ocp_qp_dims_get_ni(dims, i);
        blasfeo_dvecnrm_inf(2*ni_stage, qp_out->lam+i, 0, &res_stage);
        res = res > res_stage ? res : res_stage;
    }
    return res;
}


void ocp_qp_out_copy(ocp_qp_out* from, ocp_qp_out* to)
{
    d_ocp_qp_sol_copy_all(from, to);
}


void ocp_qp_out_axpy(double alpha, ocp_qp_out* x, ocp_qp_out* y, ocp_qp_out* z)
{
    ocp_qp_dims *dims = x->dim;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int ni_stage;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        ni_stage = ocp_qp_dims_get_ni(dims, i);
        blasfeo_daxpy(nx[i]+nu[i]+2*ns[i], alpha, x->ux+i, 0, y->ux+i, 0, z->ux+i, 0);
        blasfeo_daxpy(2*ni_stage, alpha, x->lam+i, 0, y->lam+i, 0, z->lam+i, 0);
        if (i < N)
        {
            blasfeo_daxpy(nx[i+1], alpha, x->pi+i, 0, y->pi+i, 0, z->pi+i, 0);
        }
        blasfeo_daxpy(2*ni_stage, alpha, x->t+i, 0, y->t+i, 0, z->t+i, 0);
    }
}


double ocp_qp_out_ddot(ocp_qp_out *x, ocp_qp_out *y, struct blasfeo_dvec *work_tmp_2ni)
{
    ocp_qp_dims *dims = x->dim;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nbx = dims->nbx;
    int *nbu = dims->nbu;
    int *ng = dims->ng;
    int *ns = dims->ns;
    int tmp_nbg;
    double out = 0.0;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // primal
        /* standard variant */
#if 1
        out += blasfeo_ddot(nx[i]+nu[i]+2*ns[i], x->ux+i, 0, y->ux+i, 0);
#else
        // split variant
        out += blasfeo_ddot(nx[i]+nu[i], x->ux+i, 0, y->ux+i, 0);
        // slacks
        blasfeo_daxpy(ns[i], -1.0, x->ux+i, nx[i]+nu[i], x->ux+i, nx[i]+nu[i]+ns[i], work_tmp_2ni, 0);
        blasfeo_daxpy(ns[i], -1.0, y->ux+i, nx[i]+nu[i], y->ux+i, nx[i]+nu[i]+ns[i], work_tmp_2ni, ns[i]);
        out += blasfeo_ddot(ns[i], work_tmp_2ni, 0, work_tmp_2ni, ns[i]);
#endif

        // dual
        tmp_nbg = nbu[i]+nbx[i]+ng[i];
        /* setup multipliers as lower - upper bound */
        // work_tmp_2ni = x.lam_lower - x.lam_upper
        blasfeo_daxpy(tmp_nbg, -1.0, x->lam+i, tmp_nbg, x->lam+i, 0, work_tmp_2ni, 0);
        // work_tmp_2ni[nbg:] = y.lam_lower - y.lam_upper
        blasfeo_daxpy(tmp_nbg, -1.0, y->lam+i, tmp_nbg, y->lam+i, 0, work_tmp_2ni, tmp_nbg);
        /* compute ddot of split multipliers*/
        // add dot product
        out += blasfeo_ddot(tmp_nbg, work_tmp_2ni, 0, work_tmp_2ni, tmp_nbg);
        // multipliers wrt slack bounds
#if 1
        out += blasfeo_ddot(2*ns[i], x->lam+i, 2*tmp_nbg, y->lam+i, 2*tmp_nbg);
#else
        // split variant
        blasfeo_daxpy(ns[i], -1.0, x->lam+i, 2*tmp_nbg, x->lam+i, 2*tmp_nbg+ns[i], work_tmp_2ni, 0);
        blasfeo_daxpy(ns[i], -1.0, y->lam+i, 2*tmp_nbg, y->lam+i, 2*tmp_nbg+ns[i], work_tmp_2ni, ns[i]);
        out += blasfeo_ddot(ns[i], work_tmp_2ni, 0, work_tmp_2ni, ns[i]);
#endif
        if (i < N)
        {
            out += blasfeo_ddot(nx[i+1], x->pi+i, 0, y->pi+i, 0);
        }
    }
    return out;
}


void ocp_qp_out_set_to_zero(ocp_qp_out* x)
{
    ocp_qp_dims *dims = x->dim;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int ni_stage;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        ni_stage = ocp_qp_dims_get_ni(dims, i);
        blasfeo_dvecse(nx[i]+nu[i]+2*ns[i], 0.0, x->ux+i, 0);
        blasfeo_dvecse(2*ni_stage, 0.0, x->lam+i, 0);
        if (i < N)
        {
            blasfeo_dvecse(nx[i+1], 0.0, x->pi+i, 0);
        }
        blasfeo_dvecse(2*ni_stage, 0.0, x->t+i, 0);
    }
}


void ocp_qp_out_sc(double alpha, ocp_qp_out* x)
{
    ocp_qp_dims *dims = x->dim;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int ni_stage;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        ni_stage = ocp_qp_dims_get_ni(dims, i);
        blasfeo_dvecsc(nx[i]+nu[i]+2*ns[i], alpha, x->ux+i, 0);
        blasfeo_dvecsc(2*ni_stage, alpha, x->lam+i, 0);
        if (i < N)
        {
            blasfeo_dvecsc(nx[i+1], alpha, x->pi+i, 0);
        }
        blasfeo_dvecsc(2*ni_stage, alpha, x->t+i, 0);
    }
}


// y += alpha * x
void ocp_qp_out_add(double alpha, ocp_qp_out* x, ocp_qp_out* y)
{
    ocp_qp_dims *dims = x->dim;
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    int ni_stage;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        ni_stage = ocp_qp_dims_get_ni(dims, i);
        blasfeo_dvecad(nx[i]+nu[i]+2*ns[i], alpha, x->ux+i, 0, y->ux+i, 0);
        blasfeo_dvecad(2*ni_stage, alpha, x->lam+i, 0, y->lam+i, 0);
        if (i < N)
        {
            blasfeo_dvecad(nx[i+1], alpha, x->pi+i, 0, y->pi+i, 0);
        }
        blasfeo_dvecad(2*ni_stage, alpha, x->t+i, 0, y->t+i, 0);
    }
}



/************************************************
 * res
 ************************************************/

acados_size_t ocp_qp_res_calculate_size(ocp_qp_dims *dims)
{
    acados_size_t size = sizeof(ocp_qp_res);
    size += d_ocp_qp_res_memsize(dims);
    return size;
}



ocp_qp_res *ocp_qp_res_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_res *qp_res = (ocp_qp_res *) c_ptr;
    c_ptr += sizeof(ocp_qp_res);

    d_ocp_qp_res_create(dims, qp_res, c_ptr);
    c_ptr += d_ocp_qp_res_memsize(dims);

    assert((char *) raw_memory + ocp_qp_res_calculate_size(dims) == c_ptr);

    return qp_res;
}



acados_size_t ocp_qp_res_workspace_calculate_size(ocp_qp_dims *dims)
{
    acados_size_t size = sizeof(ocp_qp_res_ws);
    size += d_ocp_qp_res_ws_memsize(dims);
    return size;
}



ocp_qp_res_ws *ocp_qp_res_workspace_assign(ocp_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_res_ws *qp_res_ws = (ocp_qp_res_ws *) c_ptr;
    c_ptr += sizeof(ocp_qp_res_ws);

    d_ocp_qp_res_ws_create(dims, qp_res_ws, c_ptr);
    c_ptr += d_ocp_qp_res_ws_memsize(dims);

    assert((char *) raw_memory + ocp_qp_res_workspace_calculate_size(dims) == c_ptr);

    return qp_res_ws;
}



void ocp_qp_res_compute(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_res *qp_res,
                        ocp_qp_res_ws *res_ws)
{

    qp_info *info = (qp_info *) qp_out->misc;

    if (info->t_computed == 0)
    {
        ocp_qp_compute_t(qp_in, qp_out);
        info->t_computed = 1;
    }

    d_ocp_qp_res_compute(qp_in, qp_out, qp_res, res_ws);

    /* mask out disregarded constraints; NOTE: not needed anymore as done in HPIPM */
    // int *nx = qp_res->dim->nx;
    // int *nu = qp_res->dim->nu;
    // int *nb = qp_res->dim->nb;
    // int *ng = qp_res->dim->ng;
    // int *ns = qp_res->dim->ns;
    // int ni_stage;
    // for (int ii = 0; ii <= qp_res->dim->N; ii++)
    // {
    //     ni_stage = ocp_qp_dims_get_ni(qp_res->dim, ii);
    //     // stationarity wrt slacks
    //     blasfeo_dvecmul(2*ns[ii], qp_in->d_mask+ii, 2*nb[ii]+2*ng[ii], &qp_res->res_g[ii], nx[ii]+nu[ii], &qp_res->res_g[ii], nx[ii]+nu[ii]);
    //     // ineq
    //     blasfeo_dvecmul(2*ni_stage, qp_in->d_mask+ii, 0, &qp_res->res_d[ii], 0, &qp_res->res_d[ii], 0);
    //     // comp
    //     blasfeo_dvecmul(2*ni_stage, qp_in->d_mask+ii, 0, &qp_res->res_m[ii], 0, &qp_res->res_m[ii], 0);
    // }
    return;
}



void ocp_qp_res_compute_nrm_inf(ocp_qp_res *qp_res, double res[4])
{
    // loop index
    int ii;

    // extract ocp qp size
    int N = qp_res->dim->N;
    int *nx = qp_res->dim->nx;
    int *nu = qp_res->dim->nu;
    int *nb = qp_res->dim->nb;
    int *ng = qp_res->dim->ng;
    int *ns = qp_res->dim->ns;

#if 1

    // compute infinity norms of residuals
    double tmp;

    res[0] = 0.0;
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_dvecnrm_inf(nx[ii] + nu[ii] + 2 * ns[ii], &qp_res->res_g[ii], 0, &tmp);
        res[0] = tmp > res[0] ? tmp : res[0];
    }

    res[1] = 0.0;
    for (ii = 0; ii < N; ii++)
    {
        blasfeo_dvecnrm_inf(nx[ii + 1], &qp_res->res_b[ii], 0, &tmp);
        res[1] = tmp > res[1] ? tmp : res[1];
    }

    res[2] = 0.0;
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_dvecnrm_inf(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_res->res_d[ii], 0, &tmp);
        res[2] = tmp > res[2] ? tmp : res[2];
    }

    res[3] = 0.0;
    for (ii = 0; ii <= N; ii++)
    {
        blasfeo_dvecnrm_inf(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_res->res_m[ii], 0, &tmp);
        res[3] = tmp > res[3] ? tmp : res[3];
    }

#else

    // XXX this should be avoided, since it does employ strucutre of the HPIPM core that may change
    // !!!

    // compute size of res_q, res_b, res_d and res_m
    int nvt = 0;
    int net = 0;
    int nct = 0;

    for (ii = 0; ii < N; ii++)
    {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];

    // compute infinity norms of residuals
    blasfeo_dvecnrm_inf(nvt, qp_res->res_g, 0, &res[0]);
    blasfeo_dvecnrm_inf(net, qp_res->res_b, 0, &res[1]);
    blasfeo_dvecnrm_inf(nct, qp_res->res_d, 0, &res[2]);
    blasfeo_dvecnrm_inf(nct, qp_res->res_m, 0, &res[3]);

#endif

    return;
}



void ocp_qp_stack_slacks_dims(ocp_qp_dims *in, ocp_qp_dims *out)
{
    int N     = in->N;
    int *nx   = in->nx;
    int *nu   = in->nu;
    int *nb   = in->nb;
    int *nbx  = in->nbx;
    int *nbu  = in->nbu;
    int *ng   = in->ng;
    int *ns   = in->ns;
    int *nsbx = in->nsbx;
    int *nsbu = in->nsbu;
    // int *nsg  = in->nsg;

    out->N = N;

    for (int ii = 0; ii <= N; ii++)
    {
        out->nx[ii] = nx[ii];
        out->nu[ii] = nu[ii] + 2*ns[ii];
        out->nb[ii] = nb[ii] - nsbx[ii] - nsbu[ii] + 2*ns[ii];
        out->nbx[ii] = nbx[ii] - nsbx[ii];
        out->nbu[ii] = nbu[ii] - nsbu[ii] + 2*ns[ii];
        out->ng[ii] = ns[ii] > 0 ? ng[ii] + nsbx[ii] + nsbu[ii] : ng[ii];
        out->ns[ii] = 0;
        out->nsbx[ii] = 0;
        out->nsbu[ii] = 0;
        out->nsg[ii] = 0;
    }
}



void ocp_qp_stack_slacks(ocp_qp_in *in, ocp_qp_in *out)
{
    int N     = in->dim->N;
    int *nx   = in->dim->nx;
    int *nu   = in->dim->nu;
    int *nb   = in->dim->nb;
    // int *nbx  = in->dim->nbx;
    // int *nbu  = in->dim->nbu;
    int *ng   = in->dim->ng;
    int *ns   = in->dim->ns;
    int *nsbx = in->dim->nsbx;
    int *nsbu = in->dim->nsbu;
    // int *nsg  = in->dim->nsg;
    int **idxb = in->idxb;
    int **idxs_rev = in->idxs_rev;

    int *nx2  = out->dim->nx;
    int *nu2  = out->dim->nu;
    int *nb2  = out->dim->nb;
    // int *nbx2 = out->dim->nbx;
    // int *nbu2 = out->dim->nbu;
    int *ng2  = out->dim->ng;

    UNUSED(nsbx);
    UNUSED(nsbu);

    for (int ii = 0; ii <= N; ii++)
    {
        if (ii < N)
        {
            // set matrices to 0.0
            blasfeo_dgese(nu2[ii]+nx2[ii]+1, nx2[ii+1], 0.0, out->BAbt+ii, 0, 0);

            // copy in->BAbt to out->BAbt, out->BAbt = [0; B'; A'; b']
            blasfeo_dgecp(nu[ii]+nx[ii]+1, nx[ii+1], in->BAbt+ii, 0, 0, out->BAbt+ii, 2*ns[ii], 0);

            // copy in->b to out->b
            blasfeo_dveccp(nx2[ii+1], in->b+ii, 0, out->b+ii, 0);
        }

        // set matrices to 0.0
        blasfeo_dgese(nu2[ii]+nx2[ii]+1, nu2[ii]+nx2[ii], 0.0, out->RSQrq+ii, 0, 0);
        blasfeo_dgese(nu2[ii]+nx2[ii], ng2[ii], 0.0, out->DCt+ii, 0, 0);

        // copy in->Z to the main diagonal of out->RSQrq, out->RSQrq = [Z 0 0; 0 0 0; 0 0 0; 0 0 0]
        blasfeo_ddiain(2*ns[ii], 1.0, in->Z+ii, 0, out->RSQrq+ii, 0, 0);

        // copy in->RSQrq to out->RSQrq, out->RSQrq = [Z 0 0; 0 R S; 0 S' Q; 0 r' q']
        blasfeo_dgecp(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], in->RSQrq+ii, 0, 0, out->RSQrq+ii, 2*ns[ii],
            2*ns[ii]);

        // copy in->rqz to out->rqz, out->rqz = [z r q]
        blasfeo_dveccp(2*ns[ii], in->rqz+ii, nu[ii]+nx[ii], out->rqz+ii, 0);
        blasfeo_dveccp(nu[ii]+nx[ii], in->rqz+ii, 0, out->rqz+ii, 2*ns[ii]);

        // copy in->DCt to out->DCt, out->DCt = [0 0; DCt 0]
        blasfeo_dgecp(nu[ii]+nx[ii], ng[ii], in->DCt+ii, 0, 0, out->DCt+ii, 2*ns[ii], 0);

        if (ns[ii] > 0)
        {
            // set flags for non-softened box constraints
            // use out->m temporarily for this
            for (int jj = 0; jj < nb[ii]; jj++)
            {
                // TODO(bnovoselnik): pick up some workspace for this
                BLASFEO_DVECEL(out->m+ii, jj) = 1.0;
            }

            int col_b = ng[ii];
//            for (int jj = 0; jj < ns[ii]; jj++)
//            {
//                int js = idxs[ii][jj];
            for(int js=0; js<nb[ii]+ng[ii]; js++)
            {
                int jj = idxs_rev[ii][js];
                if(jj!=-1)
                {

                    int idx_v_ls1 = jj;
                    int idx_v_us1 = ns[ii]+jj;

                    int idx_d_ls0 = js;
                    int idx_d_us0 = nb[ii]+ng[ii]+js;
                    int idx_d_ls1;
                    int idx_d_us1;

                    if (js < nb[ii])  // soft box constraint
                    {
                        // index of a soft box constraint
                        int jv = idxb[ii][js]+2*ns[ii];

                        idx_d_ls1 = nb2[ii]+col_b;
                        idx_d_us1 = 2*nb2[ii]+ng2[ii]+col_b;

                        // soft box constraint, set its flag to -1
                        BLASFEO_DVECEL(out->m+ii, js) = -1.0;

                        // insert soft box constraint into out->DCt, lb <= ux + sl - su <= ub
                        BLASFEO_DMATEL(out->DCt+ii, jv, col_b) = 1.0;
                        BLASFEO_DMATEL(out->DCt+ii, idx_v_ls1, col_b) = +1.0;
                        BLASFEO_DMATEL(out->DCt+ii, idx_v_us1, col_b) = -1.0;
                        BLASFEO_DVECEL(out->d+ii, idx_d_ls1) = BLASFEO_DVECEL(in->d+ii, idx_d_ls0);
                        BLASFEO_DVECEL(out->d+ii, idx_d_us1) = -BLASFEO_DVECEL(in->d+ii, idx_d_us0);

                        col_b++;
                    }
                    else  // soft general constraint
                    {
                        // index of a soft general constraint
                        int col_g = js - nb[ii];

                        // soft general constraint, lg <= D u + C x + sl - su <= ug
                        BLASFEO_DMATEL(out->DCt+ii, idx_v_ls1, col_g) = +1.0;
                        BLASFEO_DMATEL(out->DCt+ii, idx_v_us1, col_g) = -1.0;
                    }

                    // slacks have box constraints
                    out->idxb[ii][idx_v_ls1] = idx_v_ls1;
                    out->idxb[ii][idx_v_us1] = idx_v_us1;
                }
            }

            int k_nsb = 0;
            for (int jj = 0; jj < nb[ii]; jj++)
            {
                if (BLASFEO_DVECEL(out->m+ii, jj) > 0)
                {
                    // copy nonsoftened box constraint bounds to out->d
                    BLASFEO_DVECEL(out->d+ii, 2*ns[ii]+k_nsb) = BLASFEO_DVECEL(in->d+ii, jj);
                    BLASFEO_DVECEL(out->d+ii, nb2[ii]+ng2[ii]+2*ns[ii]+k_nsb) =
                                                        -BLASFEO_DVECEL(in->d+ii, nb[ii]+ng[ii]+jj);
                    out->idxb[ii][2*ns[ii]+k_nsb] = 2*ns[ii] + jj;
                    k_nsb++;
                }
            }

            assert(k_nsb == nb[ii] - nsbu[ii] - nsbx[ii] && "Dimensions are wrong!");

            // copy ls and us to out->lb
            blasfeo_dveccp(2*ns[ii], in->d+ii, 2*nb[ii]+2*ng[ii], out->d+ii, 0);

            // for slacks out->ub is +INFTY
            blasfeo_dvecse(2*ns[ii], 1.0e6, out->d+ii, nb2[ii]+ng2[ii]);

            // copy in->lg to out->lg
            blasfeo_dveccp(ng[ii], in->d+ii, nb[ii], out->d+ii, nb2[ii]);

            // copy in->ug to out->ug
            blasfeo_dveccpsc(ng[ii], -1.0, in->d+ii, 2*nb[ii]+ng[ii], out->d+ii, 2*nb2[ii]+ng2[ii]);

            // flip signs for ub and ug
            blasfeo_dvecsc(nb2[ii]+ng2[ii], -1.0, out->d+ii, nb2[ii]+ng2[ii]);

            // set out->m to 0.0
            blasfeo_dvecse(2*nb2[ii]+2*ng2[ii], 0.0, out->m+ii, 0);
        }
        else
        {
            blasfeo_dveccp(2*nb[ii]+2*ng[ii], in->d+ii, 0, out->d+ii, 0);
            blasfeo_dveccp(2*nb[ii]+2*ng[ii], in->m+ii, 0, out->m+ii, 0);
            for (int jj = 0; jj < nb[ii]; jj++) out->idxb[ii][jj] = in->idxb[ii][jj];
        }
    }
}



void ocp_qp_compute_t(ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    // loop index
    int ii;

    //
    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *ns = qp_in->dim->ns;

    struct blasfeo_dmat *DCt = qp_in->DCt;
    struct blasfeo_dvec *d = qp_in->d;
    int **idxb = qp_in->idxb;

    struct blasfeo_dvec *ux = qp_out->ux;
    struct blasfeo_dvec *t = qp_out->t;

    for (ii = 0; ii <= N; ii++)
    {
        // compute t slacks for bounds
        blasfeo_dvecex_sp(nb[ii], 1.0, idxb[ii], ux + ii, 0, t+ii, nb[ii] + ng[ii]);
        blasfeo_daxpby(nb[ii], 1.0, t + ii, nb[ii] + ng[ii], -1.0, d + ii, 0, t + ii, 0);
        blasfeo_daxpby(nb[ii], -1.0, t + ii, nb[ii] + ng[ii], -1.0, d + ii, nb[ii] + ng[ii], t + ii,
            nb[ii] + ng[ii]);

        // compute t slacks for general constraints
        blasfeo_dgemv_t(nu[ii] + nx[ii], ng[ii], 1.0, DCt + ii, 0, 0, ux + ii, 0, -1.0, d + ii, nb[ii],
                        t + ii, nb[ii]);
        blasfeo_dgemv_t(nu[ii] + nx[ii], ng[ii], -1.0, DCt + ii, 0, 0, ux + ii, 0, -1.0, d + ii,
                        2 * nb[ii] + ng[ii], t + ii, 2 * nb[ii] + ng[ii]);

        // compute t slacks for bounds on slack variables constraints
        for(int jj=0; jj<nb[ii]+ng[ii]; jj++)
        {
            int idx = qp_in->idxs_rev[ii][jj];
            if(idx!=-1)
            {
                BLASFEO_DVECEL(t+ii, jj) += BLASFEO_DVECEL(ux+ii, nu[ii]+nx[ii]+idx);
                BLASFEO_DVECEL(t+ii, nb[ii]+ng[ii]+jj) += BLASFEO_DVECEL(ux+ii, nu[ii]+nx[ii]+ns[ii]+idx);
            }
        }

        blasfeo_daxpy(2*ns[ii], -1.0, qp_in->d+ii, 2*nb[ii]+2*ng[ii], qp_out->ux+ii, nu[ii]+nx[ii], qp_out->t+ii, 2*nb[ii]+2*ng[ii]);
    }
}

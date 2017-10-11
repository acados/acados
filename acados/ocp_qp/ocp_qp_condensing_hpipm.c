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

#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"

#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"



int ocp_qp_condensing_hpipm_calculate_args_size(const ocp_qp_in *qp_in) {
    int N = qp_in->N;
    int size = 0;
    size += sizeof(ocp_qp_condensing_hpipm_args);
    size += (N+1)*sizeof(int);
    return size;
}



char *ocp_qp_condensing_hpipm_assign_args(const ocp_qp_in *qp_in,
    ocp_qp_condensing_hpipm_args **args, void *mem) {

    int N = qp_in->N;

    char *c_ptr = (char *) mem;

    *args = (ocp_qp_condensing_hpipm_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_hpipm_args);

    (*args)->scrapspace = c_ptr;
    c_ptr += (N+1)*sizeof(int);

    return c_ptr;
    }



static void ocp_qp_condensing_hpipm_initialize_default_args(const ocp_qp_in *qp_in,
    ocp_qp_condensing_hpipm_args *args) {

    args->res_g_max = 1e-6;
    args->res_b_max = 1e-8;
    args->res_d_max = 1e-8;
    args->res_m_max = 1e-8;
    args->iter_max = 50;
    args->alpha_min = 1e-8;
    args->mu0 = 1;

    int N = qp_in->N;

    int *ns = (int *) args->scrapspace;
    int ii;
    for (ii=0; ii < N+1; ii++) ns[ii] = 0;
    }



ocp_qp_condensing_hpipm_args *ocp_qp_condensing_hpipm_create_arguments(const ocp_qp_in *qp_in) {
    void *mem = malloc(ocp_qp_condensing_hpipm_calculate_args_size(qp_in));
    ocp_qp_condensing_hpipm_args *args;
    ocp_qp_condensing_hpipm_assign_args(qp_in, &args, mem);
    ocp_qp_condensing_hpipm_initialize_default_args(qp_in, args);

    return args;
}



int_t ocp_qp_condensing_hpipm_calculate_workspace_size(const ocp_qp_in *qp_in, void *args_) {
    return 0;
}



int_t ocp_qp_condensing_hpipm_calculate_memory_size(const ocp_qp_in *qp_in, void *args_) {
    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;
    // extract ocp qp in size
    int_t N = qp_in->N;
    int_t *nx = (int_t *)qp_in->nx;
    int_t *nu = (int_t *)qp_in->nu;
    int_t *nb = (int_t *)qp_in->nb;
    int_t *ng = (int_t *)qp_in->nc;
    int_t **hidxb = (int_t **)qp_in->idxb;

    // extract ns from args
    int_t *ns = (int_t *) args->scrapspace;

    // dummy ocp qp
    struct d_ocp_qp qp;
    qp.N = N;
    qp.nx = nx;
    qp.nu = nu;
    qp.nb = nb;
    qp.ng = ng;
    qp.idxb = hidxb;
    qp.ns = ns;

    // compute dense qp size
    int_t nvd = 0;
    int_t ned = 0;
    int_t nbd = 0;
    int_t ngd = 0;
    int_t nsd = 0;
#if 0
    // [u; x] order
    d_compute_qp_size_ocp2dense(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#else
    // [x; u] order  // XXX update with ns !!!!!
    d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#endif


    // dummy dense qp
    struct d_dense_qp qpd;
    qpd.nv = nvd;
    qpd.ne = ned;
    qpd.nb = nbd;
    qpd.ng = ngd;
    qpd.ns = nsd;

    // dummy ipm arg
    struct d_dense_qp_ipm_arg ipm_arg;
    ipm_arg.stat_max = args->iter_max;

    // size in bytes
    int_t size = sizeof(ocp_qp_condensing_hpipm_memory);

    size += 1 * sizeof(struct d_ocp_qp);                       // qp
    size += 1 * sizeof(struct d_ocp_qp_sol);                   // qp_sol
    size += 1 * sizeof(struct d_dense_qp);                     // qpd
    size += 1 * sizeof(struct d_dense_qp_sol);                 // qpd_sol
    size += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);  // cond_workspace
    size += 1 * sizeof(struct d_dense_qp_ipm_arg);        // ipm_arg
    size += 1 * sizeof(struct d_dense_qp_ipm_workspace);  // ipm_workspace

    size += d_memsize_ocp_qp(N, nx, nu, nb, ng, ns);
    size += d_memsize_ocp_qp_sol(N, nx, nu, nb, ng, ns);
    size += d_memsize_dense_qp(nvd, ned, nbd, ngd, nsd);
    size += d_memsize_dense_qp_sol(nvd, ned, nbd, ngd, nsd);
    size += d_memsize_cond_qp_ocp2dense(&qp, &qpd);
    size += d_memsize_dense_qp_ipm_arg(&qpd);
    size += d_memsize_dense_qp_ipm(&qpd, &ipm_arg);
    size += 4 * (N + 1) * sizeof(real_t *);  // lam_lb lam_ub lam_lg lam_ug
    size += 1 * (N + 1) * sizeof(int_t *);  // hidxb_rev
    for (int_t ii = 0; ii <= N; ii++) {
        size += nb[ii]*sizeof(int_t);  // hidxb_rev
    }

    size = (size + 63) / 64 * 64;  // make multipl of typical cache line size
    size += 1 * 64;                // align once to typical cache line size

    return size;
}



char *ocp_qp_condensing_hpipm_assign_memory(const ocp_qp_in *qp_in, void *args_, void **mem_,
                                            void *raw_memory) {

    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;
    ocp_qp_condensing_hpipm_memory **hpipm_memory = (ocp_qp_condensing_hpipm_memory **) mem_;
    // extract problem size
    int_t N = qp_in->N;
    int_t *nx = (int_t *)qp_in->nx;
    int_t *nu = (int_t *)qp_in->nu;
    int_t *nb = (int_t *)qp_in->nb;
    int_t *ng = (int_t *)qp_in->nc;
    int_t **hidxb = (int_t **)qp_in->idxb;

    // extract ns from args
    int_t *ns = (int_t *) args->scrapspace;

    // compute dense qp size
    int_t nvd = 0;
    int_t ned = 0;
    int_t nbd = 0;
    int_t ngd = 0;
    int_t nsd = 0;
#if 0
    // [u; x] order
    d_compute_qp_size_ocp2dense(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#else
    // [x; u] order  // XXX update with ns !!!!!
    d_compute_qp_size_ocp2dense_rev(N, nx, nu, nb, hidxb, ng, &nvd, &ned, &nbd, &ngd);
#endif


    // char pointer
    char *c_ptr = (char *)raw_memory;

    *hpipm_memory = (ocp_qp_condensing_hpipm_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_condensing_hpipm_memory);

    //
    (*hpipm_memory)->qp = (struct d_ocp_qp *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp);
    //
    (*hpipm_memory)->qp_sol = (struct d_ocp_qp_sol *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp_sol);
    //
    (*hpipm_memory)->qpd = (struct d_dense_qp *)c_ptr;
    c_ptr += 1 * sizeof(struct d_dense_qp);
    //
    (*hpipm_memory)->qpd_sol = (struct d_dense_qp_sol *)c_ptr;
    c_ptr += 1 * sizeof(struct d_dense_qp_sol);
    //
    (*hpipm_memory)->cond_workspace =
        (struct d_cond_qp_ocp2dense_workspace *)c_ptr;
    c_ptr += 1 * sizeof(struct d_cond_qp_ocp2dense_workspace);
    //
    (*hpipm_memory)->ipm_arg = (struct d_dense_qp_ipm_arg *)c_ptr;
    c_ptr += 1 * sizeof(struct d_dense_qp_ipm_arg);
    //
    (*hpipm_memory)->ipm_workspace = (struct d_dense_qp_ipm_workspace *)c_ptr;
    c_ptr += 1 * sizeof(struct d_dense_qp_ipm_workspace);
    //
    (*hpipm_memory)->hlam_lb = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);
    //
    (*hpipm_memory)->hlam_ub = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);
    //
    (*hpipm_memory)->hlam_lg = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);
    //
    (*hpipm_memory)->hlam_ug = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);
    //
    (*hpipm_memory)->hidxb_rev = (int_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(int_t *);

    //
    struct d_ocp_qp *qp = (*hpipm_memory)->qp;
    //
    struct d_ocp_qp_sol *qp_sol = (*hpipm_memory)->qp_sol;
    //
    struct d_dense_qp *qpd = (*hpipm_memory)->qpd;
    //
    struct d_dense_qp_sol *qpd_sol = (*hpipm_memory)->qpd_sol;
    //
    struct d_cond_qp_ocp2dense_workspace *cond_workspace =
        (*hpipm_memory)->cond_workspace;
    //
    struct d_dense_qp_ipm_arg *ipm_arg = (*hpipm_memory)->ipm_arg;
    //
    struct d_dense_qp_ipm_workspace *ipm_workspace =
        (*hpipm_memory)->ipm_workspace;

    // align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    // ocp qp structure
    d_create_ocp_qp(N, nx, nu, nb, ng, ns, qp, c_ptr);
    c_ptr += qp->memsize;
    // ocp qp sol structure
    d_create_ocp_qp_sol(N, nx, nu, nb, ng, ns, qp_sol, c_ptr);
    c_ptr += qp_sol->memsize;
    // dense qp structure
    d_create_dense_qp(nvd, ned, nbd, ngd, nsd, qpd, c_ptr);
    c_ptr += qpd->memsize;
    // dense qp sol structure
    d_create_dense_qp_sol(nvd, ned, nbd, ngd, nsd, qpd_sol, c_ptr);
    c_ptr += qpd_sol->memsize;
    // cond workspace structure
    d_create_cond_qp_ocp2dense(qp, qpd, cond_workspace, c_ptr);
    c_ptr += cond_workspace->memsize;
    // ipm arg structure
    d_create_dense_qp_ipm_arg(qpd, ipm_arg, c_ptr);
    c_ptr += ipm_arg->memsize;
    d_set_default_dense_qp_ipm_arg(ipm_arg);
    ipm_arg->iter_max = args->iter_max;
    ipm_arg->stat_max = args->iter_max;
    ipm_arg->alpha_min = args->alpha_min;
    ipm_arg->res_g_max = args->res_g_max;
    ipm_arg->res_b_max = args->res_b_max;
    ipm_arg->res_d_max = args->res_d_max;
    ipm_arg->res_m_max = args->res_m_max;
    ipm_arg->mu0 = args->mu0;
    // ipm workspace structure
    d_create_dense_qp_ipm(qpd, ipm_arg, ipm_workspace, c_ptr);
    c_ptr += ipm_workspace->memsize;

    //
    for (int_t ii = 0; ii <= N; ii++) {
        (*hpipm_memory)->hidxb_rev[ii] = (int_t *) c_ptr;
        c_ptr += nb[ii]*sizeof(int_t);
    }

    return c_ptr;
}



ocp_qp_condensing_hpipm_memory *ocp_qp_condensing_hpipm_create_memory(const ocp_qp_in *qp_in,
                                                                      void *args_) {

    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;

    ocp_qp_condensing_hpipm_memory *mem;
    int_t memory_size = ocp_qp_condensing_hpipm_calculate_memory_size(qp_in, args);
    void *raw_memory = malloc(memory_size);
    char *ptr_end =
        ocp_qp_condensing_hpipm_assign_memory(qp_in, args, (void **) &mem, raw_memory);
    assert((char*) raw_memory + memory_size >= ptr_end); (void) ptr_end;

    return mem;
}



int_t ocp_qp_condensing_hpipm(const ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                            void *args_,
                            void *memory_,
                            void *workspace_) {

    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;
    ocp_qp_condensing_hpipm_memory *memory = (ocp_qp_condensing_hpipm_memory *) memory_;

    // initialize return code
    int_t acados_status = ACADOS_SUCCESS;

    // loop index
    int_t ii, jj;

    // extract memory
    real_t **hlam_lb = memory->hlam_lb;
    real_t **hlam_ub = memory->hlam_ub;
    real_t **hlam_lg = memory->hlam_lg;
    real_t **hlam_ug = memory->hlam_ug;
    struct d_ocp_qp *qp = memory->qp;
    struct d_ocp_qp_sol *qp_sol = memory->qp_sol;
    struct d_dense_qp *qpd = memory->qpd;
    struct d_dense_qp_sol *qpd_sol = memory->qpd_sol;
    struct d_cond_qp_ocp2dense_workspace *cond_workspace =
        memory->cond_workspace;
    struct d_dense_qp_ipm_arg *ipm_arg = memory->ipm_arg;
    struct d_dense_qp_ipm_workspace *ipm_workspace = memory->ipm_workspace;
    int_t **hidxb_rev = (int_t **) memory->hidxb_rev;

    // extract problem size
    int_t N = qp_in->N;
    int_t *nx = (int_t *) qp_in->nx;
    int_t *nu = (int_t *) qp_in->nu;
    int_t *nb = (int_t *)qp_in->nb;
    int_t *ng = (int_t *)qp_in->nc;

    // extract input data
    real_t **hA = (real_t **)qp_in->A;
    real_t **hB = (real_t **)qp_in->B;
    real_t **hb = (real_t **)qp_in->b;
    real_t **hQ = (real_t **)qp_in->Q;
    real_t **hS = (real_t **)qp_in->S;
    real_t **hR = (real_t **)qp_in->R;
    real_t **hq = (real_t **)qp_in->q;
    real_t **hr = (real_t **)qp_in->r;
    real_t **hd_lb = (real_t **)qp_in->lb;
    real_t **hd_ub = (real_t **)qp_in->ub;
    real_t **hC = (real_t **)qp_in->Cx;
    real_t **hD = (real_t **)qp_in->Cu;
    real_t **hd_lg = (real_t **)qp_in->lc;
    real_t **hd_ug = (real_t **)qp_in->uc;
    int_t **hidxb = (int_t **)qp_in->idxb;

    // extract output struct members
    real_t **hx = qp_out->x;
    real_t **hu = qp_out->u;
    real_t **hpi = qp_out->pi;
    real_t **hlam = qp_out->lam;

    // compute bounds indeces in order [u; x]
    for (ii = 0; ii <= N; ii++) {
        for (jj = 0; jj < nb[ii]; jj++) {
            if (hidxb[ii][jj] < nx[ii])  {  // state constraint
                hidxb_rev[ii][jj] = hidxb[ii][jj]+nu[ii];
            } else  {  // input constraint
                hidxb_rev[ii][jj] = hidxb[ii][jj]-nx[ii];
            }
        }
    }

    //
    for (ii = 0; ii <= N; ii++) {
        hlam_lb[ii] = hlam[ii];
        hlam_ub[ii] = hlam[ii] + nb[ii];
        hlam_lg[ii] = hlam[ii] + 2 * nb[ii];
        hlam_ug[ii] = hlam[ii] + 2 * nb[ii] + ng[ii];
    }

    // ocp qp structure
    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb_rev, hd_lb, hd_ub,
                           hC, hD, hd_lg, hd_ug, NULL, NULL, NULL, NULL, NULL, qp);

    // ocp qp sol structure
    d_cond_qp_ocp2dense(qp, qpd, cond_workspace);

    // dense qp structure


    // ipm structure

    // solve ipm
    d_solve_dense_qp_ipm(qpd, qpd_sol, ipm_arg, ipm_workspace);

    // expand solution
    d_expand_sol_dense2ocp(qp, qpd_sol, qp_sol, cond_workspace);

    // extract solution
    d_cvt_ocp_qp_sol_to_colmaj(qp, qp_sol, hu, hx, NULL, NULL, hpi,
        hlam_lb, hlam_ub, hlam_lg, hlam_ug, NULL, NULL);

    // extract iteration number
    memory->iter = ipm_workspace->iter;

    // compute infinity norm of residuals
    real_t *inf_norm_res = memory->inf_norm_res;
    real_t res_tmp;
    //
    real_t *res_g;
    inf_norm_res[0] = 0;
    res_g = ipm_workspace->res_g->pa;
    for (jj = 0; jj < qpd->nv; jj++) {
        res_tmp = fabs(res_g[jj]);
        if (res_tmp > inf_norm_res[0]) {
            inf_norm_res[0] = res_tmp;
        }
    }
    real_t *res_b;
    inf_norm_res[1] = 0;
    res_b = ipm_workspace->res_b->pa;
    for (jj = 0; jj < qpd->ne; jj++) {
        res_tmp = fabs(res_b[jj]);
        if (res_tmp > inf_norm_res[1]) {
            inf_norm_res[1] = res_tmp;
        }
    }
    real_t *res_d;
    inf_norm_res[2] = 0;
    res_d = ipm_workspace->res_d->pa;
    for (jj = 0; jj < 2 * qpd->nb + 2 * qpd->ng; jj++) {
        res_tmp = fabs(res_d[jj]);
        if (res_tmp > inf_norm_res[2]) {
            inf_norm_res[2] = res_tmp;
        }
    }
    real_t *res_m;
    inf_norm_res[3] = 0;
    res_m = ipm_workspace->res_m->pa;
    for (jj = 0; jj < 2 * qpd->nb + 2 * qpd->ng; jj++) {
        res_tmp = fabs(res_m[jj]);
        if (res_tmp > inf_norm_res[3]) {
            inf_norm_res[3] = res_tmp;
        }
    }
    inf_norm_res[4] = ipm_workspace->res_mu;

    // max number of iterations
    if (ipm_workspace->iter == args->iter_max) acados_status = ACADOS_MAXITER;
    // minimum step length
    if (ipm_workspace->stat[3 + (ipm_workspace->iter - 1) * 5] <
        args->alpha_min)
        acados_status = ACADOS_MINSTEP;

    // return
    return acados_status;
    //
}



void ocp_qp_condensing_hpipm_initialize(const ocp_qp_in *qp_in, void *args_, void **mem,
                                        void **work) {
    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;

    *mem = ocp_qp_condensing_hpipm_create_memory(qp_in, args);

    int_t work_space_size = ocp_qp_condensing_hpipm_calculate_workspace_size(qp_in, args);
    *work = malloc(work_space_size);
}



void ocp_qp_condensing_hpipm_destroy(void *mem, void *work) {
    free(mem);
    free(work);
}

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

#include "acados/ocp_qp/ocp_qp_hpipm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"

#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"



int ocp_qp_hpipm_calculate_args_size(const ocp_qp_in *qp_in) {
    int N = qp_in->N;
    int size = 0;
    size += sizeof(ocp_qp_hpipm_args);
    size += (N+1)*sizeof(int);
    return size;
}



char *ocp_qp_hpipm_assign_args(const ocp_qp_in *qp_in, ocp_qp_hpipm_args **args, void *mem) {
    int N = qp_in->N;

    char *c_ptr = (char *) mem;

    *args = (ocp_qp_hpipm_args *) c_ptr;
    c_ptr += sizeof(ocp_qp_hpipm_args);

    (*args)->scrapspace = c_ptr;
    c_ptr += (N+1)*sizeof(int);

    return c_ptr;
}



static void ocp_qp_hpipm_initialize_default_args(const ocp_qp_in *qp_in, ocp_qp_hpipm_args *args) {
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



ocp_qp_hpipm_args *ocp_qp_hpipm_create_arguments(const ocp_qp_in *qp_in) {
    void *mem = malloc(ocp_qp_hpipm_calculate_args_size(qp_in));
    ocp_qp_hpipm_args *args;
    ocp_qp_hpipm_assign_args(qp_in, &args, mem);
    ocp_qp_hpipm_initialize_default_args(qp_in, args);

    return args;
}



int ocp_qp_hpipm_calculate_workspace_size(const ocp_qp_in *qp_in, ocp_qp_hpipm_args *args) {
    return 0;
}



int ocp_qp_hpipm_calculate_memory_size(const ocp_qp_in *qp_in, ocp_qp_hpipm_args *args) {
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;

    // extract ns from args
    int_t *ns = (int_t *) args->scrapspace;

    struct d_ocp_qp qp;
    qp.N = N;
    qp.nx = nx;
    qp.nu = nu;
    qp.nb = nb;
    qp.ng = ng;
    qp.ns = ns;

    // dummy ipm arg
    struct d_ocp_qp_ipm_arg ipm_arg;
    ipm_arg.stat_max = args->iter_max;

    int size = 0;

    size += sizeof(ocp_qp_hpipm_memory);

    size += 1 * sizeof(struct d_ocp_qp);                     // qp
    size += 1 * sizeof(struct d_ocp_qp_sol);                 // qp_sol
    size += 1 * sizeof(struct d_ocp_qp_ipm_arg);        // ipm_arg
    size += 1 * sizeof(struct d_ocp_qp_ipm_workspace);  // ipm_workspace

    size += d_memsize_ocp_qp(N, nx, nu, nb, ng, ns);
    size += d_memsize_ocp_qp_sol(N, nx, nu, nb, ng, ns);
    size += d_memsize_ocp_qp_ipm_arg(&qp);
    size += d_memsize_ocp_qp_ipm(&qp, &ipm_arg);
    size += 4 * (N + 1) * sizeof(double *);  // lam_lb lam_ub lam_lg lam_ug
    size += 1 * (N + 1) * sizeof(int *);  // hidxb_rev
    for (int_t ii = 0; ii <= N; ii++) {
        size += nb[ii]*sizeof(int);  // hidxb_rev
    }

    size = (size + 63) / 64 * 64;  // make multipl of typical cache line size
    size += 1 * 64;                // align once to typical cache line size

    return size;
}



char *ocp_qp_hpipm_assign_memory(const ocp_qp_in *qp_in, ocp_qp_hpipm_args *args, void **mem_,
                                void *raw_memory) {

    ocp_qp_hpipm_memory **hpipm_memory = (ocp_qp_hpipm_memory **) mem_;

    // extract problem size
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;

    // extract ns from args
    int_t *ns = (int_t *) args->scrapspace;

    // char pointer
    char *c_ptr = (char *)raw_memory;

    *hpipm_memory = (ocp_qp_hpipm_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_hpipm_memory);

    //
    (*hpipm_memory)->qp = (struct d_ocp_qp *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp);
    //
    (*hpipm_memory)->qp_sol = (struct d_ocp_qp_sol *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp_sol);
    //
    (*hpipm_memory)->ipm_arg = (struct d_ocp_qp_ipm_arg *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp_ipm_arg);
    //
    (*hpipm_memory)->ipm_workspace = (struct d_ocp_qp_ipm_workspace *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp_ipm_workspace);
    //
    (*hpipm_memory)->hlam_lb = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*hpipm_memory)->hlam_ub = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*hpipm_memory)->hlam_lg = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*hpipm_memory)->hlam_ug = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*hpipm_memory)->hidxb_rev = (int **)c_ptr;
    c_ptr += (N + 1) * sizeof(int *);

    //
    struct d_ocp_qp *qp = (*hpipm_memory)->qp;
    //
    struct d_ocp_qp_sol *qp_sol = (*hpipm_memory)->qp_sol;
    //
    struct d_ocp_qp_ipm_arg *ipm_arg = (*hpipm_memory)->ipm_arg;
    //
    struct d_ocp_qp_ipm_workspace *ipm_workspace =
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
    // ipm arg structure
    d_create_ocp_qp_ipm_arg(qp, ipm_arg, c_ptr);
    c_ptr += ipm_arg->memsize;
    d_set_default_ocp_qp_ipm_arg(ipm_arg);
    ipm_arg->iter_max = args->iter_max;
    ipm_arg->stat_max = args->iter_max;
    ipm_arg->alpha_min = args->alpha_min;
    ipm_arg->res_g_max = args->res_g_max;
    ipm_arg->res_b_max = args->res_b_max;
    ipm_arg->res_d_max = args->res_d_max;
    ipm_arg->res_m_max = args->res_m_max;
    ipm_arg->mu0 = args->mu0;
    // ipm workspace structure
    d_create_ocp_qp_ipm(qp, ipm_arg, ipm_workspace, c_ptr);
    c_ptr += ipm_workspace->memsize;

    //
    for (int_t ii = 0; ii <= N; ii++) {
        (*hpipm_memory)->hidxb_rev[ii] = (int *) c_ptr;
        c_ptr += nb[ii]*sizeof(int);
    }

    return c_ptr;
}



ocp_qp_hpipm_memory *ocp_qp_hpipm_create_memory(const ocp_qp_in *qp_in, void *args_) {
    ocp_qp_hpipm_args *args = (ocp_qp_hpipm_args *) args_;

    ocp_qp_hpipm_memory *mem;
    int_t memory_size = ocp_qp_hpipm_calculate_memory_size(qp_in, args);
    void *raw_memory = malloc(memory_size);
    char *ptr_end = ocp_qp_hpipm_assign_memory(qp_in, args, (void **) &mem, raw_memory);
    assert((char *) raw_memory + memory_size >= ptr_end); (void) ptr_end;

    return mem;
}



int ocp_qp_hpipm(const ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args_, void *mem_, void *work_) {

    ocp_qp_hpipm_args *args = (ocp_qp_hpipm_args *) args_;
    ocp_qp_hpipm_memory *memory = (ocp_qp_hpipm_memory *) mem_;
    //
    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // loop index
    int ii, jj;

    // extract memory
    double **hlam_lb = memory->hlam_lb;
    double **hlam_ub = memory->hlam_ub;
    double **hlam_lg = memory->hlam_lg;
    double **hlam_ug = memory->hlam_ug;
    struct d_ocp_qp *qp = memory->qp;
    struct d_ocp_qp_sol *qp_sol = memory->qp_sol;
    struct d_ocp_qp_ipm_arg *ipm_arg = memory->ipm_arg;
    struct d_ocp_qp_ipm_workspace *ipm_workspace = memory->ipm_workspace;
    int **hidxb_rev = (int **) memory->hidxb_rev;

    // extract problem size
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;

    // extract input data
    double **hA = (double **)qp_in->A;
    double **hB = (double **)qp_in->B;
    double **hb = (double **)qp_in->b;
    double **hQ = (double **)qp_in->Q;
    double **hS = (double **)qp_in->S;
    double **hR = (double **)qp_in->R;
    double **hq = (double **)qp_in->q;
    double **hr = (double **)qp_in->r;
    double **hd_lb = (double **)qp_in->lb;
    double **hd_ub = (double **)qp_in->ub;
    double **hC = (double **)qp_in->Cx;
    double **hD = (double **)qp_in->Cu;
    double **hd_lg = (double **)qp_in->lc;
    double **hd_ug = (double **)qp_in->uc;
    int **hidxb = (int **)qp_in->idxb;

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam = qp_out->lam;

    // compute bounds indeces in order [u; x]
    for (ii = 0; ii <= N; ii++) {
        for (jj = 0; jj < nb[ii]; jj++) {
            if (hidxb[ii][jj] < nx[ii]) {  // state constraint
                hidxb_rev[ii][jj] = hidxb[ii][jj]+nu[ii];
            } else {  // input constraint
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


    // ipm structure

    // solve ipm
    d_solve_ocp_qp_ipm(qp, qp_sol, ipm_arg, ipm_workspace);

    // extract solution
    d_cvt_ocp_qp_sol_to_colmaj(qp, qp_sol, hu, hx, NULL, NULL, hpi,
                               hlam_lb, hlam_ub, hlam_lg, hlam_ug, NULL, NULL);

    // extract iteration number
    memory->iter = ipm_workspace->iter;

    // compute infinity norm of residuals
    double *inf_norm_res = memory->inf_norm_res;
    double res_tmp;
    //
    double *res_g;
    inf_norm_res[0] = 0;
    for (ii = 0; ii <= N; ii++) {
        res_g = (ipm_workspace->res_workspace->res_g + ii)->pa;
        for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
            res_tmp = fabs(res_g[jj]);
            if (res_tmp > inf_norm_res[0]) {
                inf_norm_res[0] = res_tmp;
            }
        }
    }
    double *res_b;
    inf_norm_res[1] = 0;
    for (ii = 0; ii < N; ii++) {
        res_b = (ipm_workspace->res_workspace->res_b + ii)->pa;
        for (jj = 0; jj < nx[ii + 1]; jj++) {
            res_tmp = fabs(res_b[jj]);
            if (res_tmp > inf_norm_res[1]) {
                inf_norm_res[1] = res_tmp;
            }
        }
    }
    double *res_d;
    inf_norm_res[2] = 0;
    for (ii = 0; ii <= N; ii++) {
        res_d = (ipm_workspace->res_workspace->res_d + ii)->pa+0;
        for (jj = 0; jj < nb[ii]; jj++) {
            res_tmp = fabs(res_d[jj]);
            if (res_tmp > inf_norm_res[2]) {
                inf_norm_res[2] = res_tmp;
            }
        }
        res_d = (ipm_workspace->res_workspace->res_d + ii)->pa+nb[ii]+ng[ii];
        for (jj = 0; jj < nb[ii]; jj++) {
            res_tmp = fabs(res_d[jj]);
            if (res_tmp > inf_norm_res[2]) {
                inf_norm_res[2] = res_tmp;
            }
        }
        res_d = (ipm_workspace->res_workspace->res_d + ii)->pa+nb[ii];
        for (jj = 0; jj < ng[ii]; jj++) {
            res_tmp = fabs(res_d[jj]);
            if (res_tmp > inf_norm_res[2]) {
                inf_norm_res[2] = res_tmp;
            }
        }
        res_d = (ipm_workspace->res_workspace->res_d + ii)->pa+2*nb[ii]+ng[ii];
        for (jj = 0; jj < ng[ii]; jj++) {
            res_tmp = fabs(res_d[jj]);
            if (res_tmp > inf_norm_res[2]) {
                inf_norm_res[2] = res_tmp;
            }
        }
    }
    double *res_m;
    inf_norm_res[3] = 0;
    for (ii = 0; ii <= N; ii++) {
        res_m = (ipm_workspace->res_workspace->res_m + ii)->pa+0;
        for (jj = 0; jj < nb[ii]; jj++) {
            res_tmp = fabs(res_m[jj]);
            if (res_tmp > inf_norm_res[3]) {
                inf_norm_res[3] = res_tmp;
            }
        }
        res_m = (ipm_workspace->res_workspace->res_m + ii)->pa+nb[ii]+ng[ii];
        for (jj = 0; jj < nb[ii]; jj++) {
            res_tmp = fabs(res_m[jj]);
            if (res_tmp > inf_norm_res[3]) {
                inf_norm_res[3] = res_tmp;
            }
        }
        res_m = (ipm_workspace->res_workspace->res_m + ii)->pa+nb[ii];
        for (jj = 0; jj < ng[ii]; jj++) {
            res_tmp = fabs(res_m[jj]);
            if (res_tmp > inf_norm_res[3]) {
                inf_norm_res[3] = res_tmp;
            }
        }
        res_m = (ipm_workspace->res_workspace->res_m + ii)->pa+2*nb[ii]+ng[ii];
        for (jj = 0; jj < ng[ii]; jj++) {
            res_tmp = fabs(res_m[jj]);
            if (res_tmp > inf_norm_res[3]) {
                inf_norm_res[3] = res_tmp;
            }
        }
    }
    inf_norm_res[4] = ipm_workspace->res_workspace->res_mu;

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



void ocp_qp_hpipm_initialize(const ocp_qp_in *qp_in, void *args_, void **mem, void **work) {
    ocp_qp_hpipm_args *args = (ocp_qp_hpipm_args *) args_;

    *mem = ocp_qp_hpipm_create_memory(qp_in, args);

    int_t work_space_size = ocp_qp_hpipm_calculate_workspace_size(qp_in, args);
    *work = malloc(work_space_size);
}

void ocp_qp_hpipm_destroy(void *mem, void *work) {
    free(mem);
    free(work);
}

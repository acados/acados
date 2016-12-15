/*
 *    This file is part of ACADOS.
 *
 *    ACADOS is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADOS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADOS; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <stdlib.h>

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#include "hpmpc/include/c_interface.h"

#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_blas.h>
#include <blasfeo/include/blasfeo_d_aux.h>

// work space size
int ocp_qp_hpmpc_workspace_size_bytes(int N, int *nx, int *nu, int *nb, int *ng, int **hidxb, \
    ocp_qp_hpmpc_args *hpmpc_args) {

    int ii;

    int N2 = hpmpc_args->N2;

    int k_max = hpmpc_args->max_iter;

    int workspace_size = 8 + 5*k_max*sizeof(double);  // alignment to single-word (4-byte)

    for (ii=0; ii <= N; ii++)
        workspace_size += nb[ii]*sizeof(int);

    workspace_size += hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(N, nx, nu, nb, hidxb, ng, N2);

    return workspace_size;
}

int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_hpmpc_args *hpmpc_args, \
    void *workspace) {

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // loop index
    int ii, jj;

    // extract input struct members
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int **hidxb_swp = (int **) qp_in->idxb;
    int *ng = (int *) qp_in->nc;
    double **hA = (double **) qp_in->A;
    double **hB = (double **) qp_in->B;
    double **hb = (double **) qp_in->b;
    double **hQ = (double **) qp_in->Q;
    double **hS = (double **) qp_in->S;
    double **hR = (double **) qp_in->R;
    double **hq = (double **) qp_in->q;
    double **hr = (double **) qp_in->r;
    double **hlb = (double **) qp_in->lb;
    double **hub = (double **) qp_in->ub;
    double **hC = (double **) qp_in->Cx;
    double **hD = (double **) qp_in->Cu;
    double **hlg = (double **) qp_in->lc;
    double **hug = (double **) qp_in->uc;

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam = qp_out->lam;

    // extract args struct members
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;
    int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem

    //  other solver arguments
    int kk = -1;  // actual number of iterations
    double inf_norm_res[4];  // inf norm of residuals
    for (ii = 0; ii < 4; ii++) inf_norm_res[ii] = 0.0;  // zero

    // memory for stat
    size_t addr = ( ( (size_t) workspace) + 7) / 8 * 8;  // align to 8-byte boundaries
    double *ptr_double = (double *) addr;
    double *stat = ptr_double;
    ptr_double += 5*k_max;
    for (ii = 0; ii < 5*k_max; ii++) stat[ii] = 0.0;  // zero

    // memory for idxb
    int *hidxb[N+1];
//  size_t addr = (( (size_t) workspace ) + 3 ) / 4 * 4;  // align to 4-byte boundaries
    int *ptr_int = (int *) ptr_double;
    for (ii = 0; ii <= N; ii++) {
        hidxb[ii] = ptr_int;
        ptr_int += nb[ii];
    }
    workspace = (void *) ptr_int;

    //  swap x and u in bounds (by updating their indeces)
    for (ii = 0; ii <= N; ii++) {
        jj = 0;
        for (; jj < nb[ii]; jj++) {
            if (hidxb_swp[ii][jj] < nx[ii]) {  // state
                hidxb[ii][jj] = hidxb_swp[ii][jj]+nu[ii];
            } else {  // input
                hidxb[ii][jj] = hidxb_swp[ii][jj]-nx[ii];
            }
        }
    }

    int hpmpc_status = fortran_order_d_ip_ocp_hard_tv(&kk, k_max, mu0, mu_tol, N, nx, nu, nb, \
        hidxb, ng, N2, warm_start, hA, hB, hb, hQ, hS, hR, hq, hr, hlb, hub, hC, hD, hlg, hug, \
        hx, hu, hpi, hlam, inf_norm_res, workspace, stat);

    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;

    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;

    // return
    return acados_status;
}

int ocp_qp_hpmpc_libstr(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_hpmpc_args *hpmpc_args, \
    void *workspace) {

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // loop index
    int ii, jj;

    // extract input struct members
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int **hsidxb = (int **) qp_in->idxb;
    int *ng = (int *) qp_in->nc;
    double **hA = (double **) qp_in->A;
    double **hB = (double **) qp_in->B;
    double **hb = (double **) qp_in->b;
    double **hQ = (double **) qp_in->Q;
    double **hS = (double **) qp_in->S;
    double **hR = (double **) qp_in->R;
    double **hq = (double **) qp_in->q;
    double **hr = (double **) qp_in->r;
    double **hlb = (double **) qp_in->lb;
    double **hub = (double **) qp_in->ub;
    double **hC = (double **) qp_in->Cx;
    double **hD = (double **) qp_in->Cu;
    double **hlg = (double **) qp_in->lc;
    double **hug = (double **) qp_in->uc;

    struct d_strmat hsBAbt[N+1];
    struct d_strvec hsb[N+1];
    struct d_strmat hsRSQrq[N+1];
    struct d_strvec hsrq[N+1];
    struct d_strmat hsDCt[N+1];
    struct d_strvec hsd[N+1];
    struct d_strvec hsux[N+1];
    struct d_strvec hspi[N+1];
    struct d_strvec hslam[N+1];
    struct d_strvec hst[N+1];
    struct d_strvec hsPb[N+1];
    struct d_strmat hsL[N+1];
    struct d_strmat hsLxt[N+1];
    struct d_strmat hsric_work_mat[2];
    struct d_strvec hsric_work_vec[1];

    int nuM;
    int nbM;

    char *ptr_memory = (char *) workspace;

    for( ii=0; ii<N; ii++ ) {
      d_create_strmat(nu[ii]+nx[ii]+1, nx[ii+1], &hsBAbt[ii+1], ptr_memory);
      ptr_memory += (&hsBAbt[ii+1])->memory_size;
  		d_cvt_tran_mat2strmat(nx[ii+1], nu[ii], hB[ii], nx[ii+1], &hsBAbt[ii+1], 0, 0);
  		d_cvt_tran_mat2strmat(nx[ii+1], nx[ii], hA[ii], nx[ii+1], &hsBAbt[ii+1], nu[ii], 0);
    	d_cvt_tran_mat2strmat(nx[ii+1], 1, hb[ii], nx[ii+1], &hsBAbt[ii+1], nu[ii]+nx[ii], 0);

      d_create_strvec(nx[ii+1], &hsb[ii+1], ptr_memory);
      ptr_memory += (&hsb[ii+1])->memory_size;
      d_cvt_vec2strvec(nx[ii+1], hb[ii], &hsb[ii+1], 0);

      d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], ptr_memory);
      ptr_memory += (&hsRSQrq[ii])->memory_size;
      d_cvt_mat2strmat(nu[ii], nu[ii], hR[ii], nu[ii], &hsRSQrq[ii], 0, 0);
      d_cvt_tran_mat2strmat(nu[ii], nx[ii], hS, nu[ii], &hsRSQrq[ii], nu[ii], 0);
      d_cvt_mat2strmat(nx[ii], nx[ii], hQ[ii], nx[ii], &hsRSQrq[ii], nu[ii], nu[ii]);
      d_cvt_tran_mat2strmat(nu[ii], 1, hr[ii], nu[ii], &hsRSQrq[ii], nu[ii]+nx[ii], 0);
      d_cvt_tran_mat2strmat(nx[ii], 1, hq[ii], nx[ii], &hsRSQrq[ii], nu[ii]+nx[ii], nu[ii]);

      d_create_strvec(nu[ii]+nx[ii], &hsrq[ii], ptr_memory);
      ptr_memory += (&hsrq[ii])->memory_size;
      d_cvt_vec2strvec(nu[ii], hr, &hsrq[ii], 0);
      d_cvt_vec2strvec(nx[ii], hq, &hsrq[ii], nu[ii]);

      d_create_strmat(nu[ii]+nx[ii]+1, ng[ii], &hsDCt[ii], ptr_memory);
      ptr_memory += (&hsDCt[ii])->memory_size;
      d_cvt_tran_mat2strmat(ng[ii], nu[ii], hD[ii], ng[ii], &hsDCt[ii], 0, 0);
      d_cvt_tran_mat2strmat(ng[ii], nx[ii], hC[ii], ng[ii], &hsDCt[ii], nu[ii], 0);

      d_create_strvec(2*nb[ii]+2*ng[ii], &hsd[ii], ptr_memory);
      ptr_memory += (&hsd[ii])->memory_size;
      d_cvt_vec2strvec(nb[ii], hlb[ii], &hsd[ii], 0);
      d_cvt_vec2strvec(nb[ii], hub[ii], &hsd[ii], nb[ii]);
      d_cvt_vec2strvec(ng[ii], hlg[ii], &hsd[ii], 2*nb[ii]);
      d_cvt_vec2strvec(ng[ii], hug[ii], &hsd[ii], 2*nb[ii] + ng[ii]);

      d_create_strvec(nu[ii]+nx[ii], &hsux[ii], ptr_memory);
      ptr_memory += (&hsux[ii])->memory_size;

      d_create_strvec(nx[ii+1], &hspi[ii+1], ptr_memory);
      ptr_memory += (&hspi[ii+1])->memory_size;

      d_create_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
      ptr_memory += (&hslam[ii])->memory_size;

      d_create_strvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
      ptr_memory += (&hst[ii])->memory_size;
    }

    ii = N;

    d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], ptr_memory);
    ptr_memory += (&hsRSQrq[ii])->memory_size;
    d_cvt_mat2strmat(nu[ii], nu[ii], hR[ii], nu[ii], &hsRSQrq[ii], 0, 0);
    d_cvt_tran_mat2strmat(nu[ii], nx[ii], hS, nu[ii], &hsRSQrq[ii], nu[ii], 0);
    d_cvt_mat2strmat(nx[ii], nx[ii], hQ[ii], nx[ii], &hsRSQrq[ii], nu[ii], nu[ii]);
    d_cvt_tran_mat2strmat(nu[ii], 1, hr[ii], nu[ii], &hsRSQrq[ii], nu[ii]+nx[ii], 0);
    d_cvt_tran_mat2strmat(nx[ii], 1, hq[ii], nx[ii], &hsRSQrq[ii], nu[ii]+nx[ii], nu[ii]);

    d_create_strvec(2*nb[ii]+2*ng[ii], &hsd[ii], ptr_memory);
    ptr_memory += (&hsd[ii])->memory_size;
    d_cvt_vec2strvec(nb[ii], hlb[ii], &hsd[ii], 0);
    d_cvt_vec2strvec(nb[ii], hub[ii], &hsd[ii], nb[ii]);
    d_cvt_vec2strvec(ng[ii], hlg[ii], &hsd[ii], 2*nb[ii]);
    d_cvt_vec2strvec(ng[ii], hug[ii], &hsd[ii], 2*nb[ii] + ng[ii]);

    d_create_strvec(nu[ii]+nx[ii], &hsux[ii], ptr_memory);
    ptr_memory += (&hsux[ii])->memory_size;

    // d_create_strvec(nx[ii+1], &hspi[ii+1], ptr_memory);
    // ptr_memory += (&hspi[ii+1])->memory_size;

    d_create_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
    ptr_memory += (&hslam[ii])->memory_size;

    d_create_strvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
    ptr_memory += (&hst[ii])->memory_size;

    // TODO: LAST STAGE?

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam = qp_out->lam;

    // extract args struct members
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;
    int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem

    //  other solver arguments
    int kk = -1;  // actual number of iterations
    double inf_norm_res[4];  // inf norm of residuals
    for (ii = 0; ii < 4; ii++) inf_norm_res[ii] = 0.0;  // zero

  	double mu;
  	struct d_strvec hswork[2];

    // max sizes
    int ngM = 0;
    for ( ii=0; ii<=N; ii++ ) {
      ngM = ng[ii]>ngM ? ng[ii] : ngM;
    }

    int nzM  = 0;
    for( ii=0; ii<=N; ii++ ) {
      nzM = nu[ii]+nx[ii]+1>nzM ? nu[ii]+nx[ii]+1 : nzM;
    }

    int nxgM = ng[N];
    for ( ii=0; ii<N; ii++ ) {
      nxgM = nx[ii+1]+ng[ii]>nxgM ? nx[ii+1]+ng[ii] : nxgM;
    }

    // IPM constants
    int hpmpc_status;
    int kk_avg;
    double alpha_min = 1e-8;
    double *stat; d_zeros(&stat, k_max, 5);
    int compute_res = 1;
    int compute_mult = 1;

    void *work_memory;
    // v_zeros_align(&work_memory,
      // d_ip2_res_mpc_hard_tv_work_space_size_bytes_libstr(N, nx, nu, nb, ng));

    hpmpc_status = d_ip2_res_mpc_hard_libstr(&kk, k_max, mu0, mu_tol, alpha_min,
      warm_start, stat, N, nx, nu, nb, hsidxb, ng, hsBAbt, hsRSQrq, hsDCt,
      hsd, hsux, compute_mult, hspi, hslam, hst, ptr_memory);

    double **temp_u;
    // copy result to qp_out
    for ( ii = 0; ii<=N; ii++ ) {
      hu[ii] = hsux[ii].pa;
      temp_u = &hsux[ii].pa;
      hx[ii] = &temp_u[0][nu[ii]];
    }

    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;

    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;

    // return
    return acados_status;
}


int ocp_qp_hpnmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_hpmpc_args *hpmpc_args, \
    void *workspace) {

    // initialize return code
    int acados_status = ACADOS_SUCCESS;

    // loop index
    int ii, jj;

    // extract input struct members
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int **hidxb_swp = (int **) qp_in->idxb;
    int *ng = (int *) qp_in->nc;
    double **hA = (double **) qp_in->A;
    double **hB = (double **) qp_in->B;
    double **hb = (double **) qp_in->b;
    double **hQ = (double **) qp_in->Q;
    double **hS = (double **) qp_in->S;
    double **hR = (double **) qp_in->R;
    double **hq = (double **) qp_in->q;
    double **hr = (double **) qp_in->r;
    double **hlb = (double **) qp_in->lb;
    double **hub = (double **) qp_in->ub;
    double **hC = (double **) qp_in->Cx;
    double **hD = (double **) qp_in->Cu;
    double **hlg = (double **) qp_in->lc;
    double **hug = (double **) qp_in->uc;

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam = qp_out->lam;
    double **ht = qp_out->t;

    // extract args struct members
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;
    int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem

    double **ux0 = hpmpc_args->ux0;
    double **pi0 = hpmpc_args->pi0;
    double **lam0 = hpmpc_args->lam0;
    double **t0 = hpmpc_args->t0;

    //  other solver arguments
    int kk = -1;  // actual number of iterations
    double inf_norm_res[4];  // inf norm of residuals
    for (ii = 0; ii < 4; ii++) inf_norm_res[ii] = 0.0;  // zero

    // memory for stat
    size_t addr = ( ( (size_t) workspace) + 7) / 8 * 8;  // align to 8-byte boundaries
    double *ptr_double = (double *) addr;
    double *stat = ptr_double;
    ptr_double += 5*k_max;
    for (ii = 0; ii < 5*k_max; ii++) stat[ii] = 0.0;  // zero

    // memory for idxb
    int *hidxb[N+1];
//  size_t addr = (( (size_t) workspace ) + 3 ) / 4 * 4;  // align to 4-byte boundaries
    int *ptr_int = (int *) ptr_double;
    for (ii = 0; ii <= N; ii++) {
        hidxb[ii] = ptr_int;
        ptr_int += nb[ii];
    }
    workspace = (void *) ptr_int;

    //  swap x and u in bounds (by updating their indeces)
    for (ii = 0; ii <= N; ii++) {
        jj = 0;
        for (; jj < nb[ii]; jj++) {
            if (hidxb_swp[ii][jj] < nx[ii]) {  // state
                hidxb[ii][jj] = hidxb_swp[ii][jj]+nu[ii];
            } else {  // input
                hidxb[ii][jj] = hidxb_swp[ii][jj]-nx[ii];
            }
        }
    }

    int hpmpc_status = fortran_order_d_ip_ocp_hard_tv_single_newton_step(&kk, k_max, mu0, mu_tol, N, nx, nu, nb, \
        hidxb, ng, N2, warm_start, hA, hB, hb, hQ, hS, hR, hq, hr, hlb, hub, hC, hD, hlg, hug, \
        hx, hu, hpi, hlam, ht, inf_norm_res, workspace, stat, ux0, pi0, lam0, t0);

    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;

    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;

    // return
    return acados_status;
}

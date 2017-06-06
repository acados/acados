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

#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "hpmpc/include/c_interface.h"
#include "hpmpc/include/mpc_solvers.h"
#include "hpmpc/include/lqcp_solvers.h"
#include "hpmpc/include/mpc_aux.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

// // work space size
// int ocp_qp_hpmpc_workspace_size_bytes(int N, int *nx, int *nu, int *nb, int *ng, int **hidxb,
//     ocp_qp_hpmpc_args *hpmpc_args) {
//
//     int ii;
//
//     int N2 = hpmpc_args->N2;
//
//     int k_max = hpmpc_args->max_iter;
//
//     int workspace_size = 8 + 5*k_max*sizeof(double);  // alignment to single-word (4-byte)
//
//     for (ii=0; ii <= N; ii++) workspace_size += nb[ii]*sizeof(int);
//         workspace_size += hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(N, nx,
//           nu, nb, hidxb, ng, N2);
//
//     return workspace_size;
// }

void ocp_qp_hpmpc_initialize(ocp_qp_in *qp_in, void *args_, void *mem_, void **work) {
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args*) args_;

    // TODO(andrea): replace dummy commands once interface completed
    args->max_iter = args->max_iter;
    if (qp_in->nx[0] > 0)
        mem_++;
    work++;
}

void ocp_qp_hpmpc_destroy(void *mem_, void *work) {
  free(work);
  ocp_qp_hpmpc_free_memory(mem_);
}

// int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
// ocp_qp_hpmpc_args *hpmpc_args,
//     void *workspace) {
//
//     // initialize return code
//     int acados_status = ACADOS_SUCCESS;
//
//     // loop index
//     int ii, jj;
//
//     // extract input struct members
//     int N = qp_in->N;
//     int *nx = (int *) qp_in->nx;
//     int *nu = (int *) qp_in->nu;
//     int *nb = (int *) qp_in->nb;
//     int **hidxb_swp = (int **) qp_in->idxb;
//     int *ng = (int *) qp_in->nc;
//     double **hA = (double **) qp_in->A;
//     double **hB = (double **) qp_in->B;
//     double **hb = (double **) qp_in->b;
//     double **hQ = (double **) qp_in->Q;
//     double **hS = (double **) qp_in->S;
//     double **hR = (double **) qp_in->R;
//     double **hq = (double **) qp_in->q;
//     double **hr = (double **) qp_in->r;
//     double **hlb = (double **) qp_in->lb;
//     double **hub = (double **) qp_in->ub;
//     double **hC = (double **) qp_in->Cx;
//     double **hD = (double **) qp_in->Cu;
//     double **hlg = (double **) qp_in->lc;
//     double **hug = (double **) qp_in->uc;
//
//     // extract output struct members
//     double **hx = qp_out->x;
//     double **hu = qp_out->u;
//     double **hpi = qp_out->pi;
//     double **hlam = qp_out->lam;
//
//     // extract args struct members
//     double mu_tol = hpmpc_args->tol;
//     int k_max = hpmpc_args->max_iter;
//     double mu0 = hpmpc_args->mu0;
//     int warm_start = hpmpc_args->warm_start;
//     int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem
//     int out_iter = -1;  // number of performed iterations
//     double *inf_norm_res = hpmpc_args->inf_norm_res;
//
//     // memory for stat
//     size_t addr = ( ( (size_t) workspace) + 7) / 8 * 8;  // align to 8-byte boundaries
//     double *ptr_double = (double *) addr;
//     double *stat = ptr_double;
//     ptr_double += 5*k_max;
//     for (ii = 0; ii < 5*k_max; ii++) stat[ii] = 0.0;  // zero
//
//     // memory for idxb
//     int *hidxb[N+1];
// //  size_t addr = (( (size_t) workspace ) + 3 ) / 4 * 4;  // align to 4-byte boundaries
//     int *ptr_int = (int *) ptr_double;
//     for (ii = 0; ii <= N; ii++) {
//         hidxb[ii] = ptr_int;
//         ptr_int += nb[ii];
//     }
//     workspace = (void *) ptr_int;
//
//     //  swap x and u in bounds (by updating their indeces)
//     for (ii = 0; ii <= N; ii++) {
//         jj = 0;
//         for (; jj < nb[ii]; jj++) {
//             if (hidxb_swp[ii][jj] < nx[ii]) {  // state
//                 hidxb[ii][jj] = hidxb_swp[ii][jj]+nu[ii];
//             } else {  // input
//                 hidxb[ii][jj] = hidxb_swp[ii][jj]-nx[ii];
//             }
//         }
//     }
//
//     int hpmpc_status = fortran_order_d_ip_ocp_hard_tv(&out_iter, k_max, mu0,
//        mu_tol, N, nx, nu, nb, hidxb, ng, N2, warm_start, hA, hB, hb, hQ, hS,
//        hR, hq, hr, hlb, hub, hC, hD, hlg,
//         hug, hx, hu, hpi, hlam, inf_norm_res, workspace, stat);
//
//     hpmpc_args->out_iter = out_iter;  // number of performed iterations
//
//     if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;
//
//     if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;
//
//     // return
//     return acados_status;
// }

// int ocp_qp_hpmpc_libstr(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
// ocp_qp_hpmpc_args *hpmpc_args, void *workspace) {
//
//     // initialize return code
//     int acados_status = ACADOS_SUCCESS;
//
//     // loop index
//     int ii;
//
//     // extract input struct members
//     int N = qp_in->N;
//     int *nx = (int *) qp_in->nx;
//     int *nu = (int *) qp_in->nu;
//     int *nb = (int *) qp_in->nb;
//     int **hsidxb = (int **) qp_in->idxb;
//     int *ng = (int *) qp_in->nc;
//     double **hA = (double **) qp_in->A;
//     double **hB = (double **) qp_in->B;
//     double **hb = (double **) qp_in->b;
//     double **hQ = (double **) qp_in->Q;
//     double **hS = (double **) qp_in->S;
//     double **hR = (double **) qp_in->R;
//     double **hq = (double **) qp_in->q;
//     double **hr = (double **) qp_in->r;
//     double **hlb = (double **) qp_in->lb;
//     double **hub = (double **) qp_in->ub;
//     double **hC = (double **) qp_in->Cx;
//     double **hD = (double **) qp_in->Cu;
//     double **hlg = (double **) qp_in->lc;
//     double **hug = (double **) qp_in->uc;
//
//     struct d_strmat hsBAbt[N+1];
//     struct d_strvec hsb[N+1];
//     struct d_strmat hsRSQrq[N+1];
//     struct d_strvec hsrq[N+1];
//     struct d_strmat hsDCt[N+1];
//     struct d_strvec hsd[N+1];
//     struct d_strvec hsux[N+1];
//     struct d_strvec hspi[N+1];
//     struct d_strvec hslam[N+1];
//     struct d_strvec hst[N+1];
//     // struct d_strvec hsPb[N+1];
//     // struct d_strmat hsL[N+1];
//     // struct d_strmat hsLxt[N+1];
//     // struct d_strmat hsric_work_mat[2];
//     // struct d_strvec hsric_work_vec[1];
//
//     // int nuM;
//     // int nbM;
//
//     char *ptr_memory = (char *) workspace;
//
//     for ( ii = 0; ii < N; ii++ ) {
//       d_create_strmat(nu[ii]+nx[ii]+1, nx[ii+1], &hsBAbt[ii], ptr_memory);
//       ptr_memory += (&hsBAbt[ii])->memory_size;
//       d_cvt_tran_mat2strmat(nx[ii+1], nu[ii], hB[ii], nx[ii+1], &hsBAbt[ii], 0, 0);
//       d_cvt_tran_mat2strmat(nx[ii+1], nx[ii], hA[ii], nx[ii+1], &hsBAbt[ii], nu[ii], 0);
//       d_cvt_tran_mat2strmat(nx[ii+1], 1, hb[ii], nx[ii+1], &hsBAbt[ii], nu[ii]+nx[ii], 0);
//
//       d_create_strvec(nx[ii+1], &hsb[ii], ptr_memory);
//       ptr_memory += (&hsb[ii])->memory_size;
//       d_cvt_vec2strvec(nx[ii+1], hb[ii], &hsb[ii], 0);
//
//       d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], ptr_memory);
//       ptr_memory += (&hsRSQrq[ii])->memory_size;
//       d_cvt_mat2strmat(nu[ii], nu[ii], hR[ii], nu[ii], &hsRSQrq[ii], 0, 0);
//       d_cvt_tran_mat2strmat(nu[ii], nx[ii], hS[ii], nu[ii], &hsRSQrq[ii], nu[ii], 0);
//       d_cvt_mat2strmat(nx[ii], nx[ii], hQ[ii], nx[ii], &hsRSQrq[ii], nu[ii], nu[ii]);
//       d_cvt_tran_mat2strmat(nu[ii], 1, hr[ii], nu[ii], &hsRSQrq[ii], nu[ii]+nx[ii], 0);
//       d_cvt_tran_mat2strmat(nx[ii], 1, hq[ii], nx[ii], &hsRSQrq[ii], nu[ii]+nx[ii], nu[ii]);
//
//       d_create_strvec(nu[ii]+nx[ii], &hsrq[ii], ptr_memory);
//       ptr_memory += (&hsrq[ii])->memory_size;
//       d_cvt_vec2strvec(nu[ii], hr[ii], &hsrq[ii], 0);
//       d_cvt_vec2strvec(nx[ii], hq[ii], &hsrq[ii], nu[ii]);
//
//       d_create_strmat(nu[ii]+nx[ii]+1, ng[ii], &hsDCt[ii], ptr_memory);
//       ptr_memory += (&hsDCt[ii])->memory_size;
//       d_cvt_tran_mat2strmat(ng[ii], nu[ii], hD[ii], ng[ii], &hsDCt[ii], 0, 0);
//       d_cvt_tran_mat2strmat(ng[ii], nx[ii], hC[ii], ng[ii], &hsDCt[ii], nu[ii], 0);
//
//       d_create_strvec(2*nb[ii]+2*ng[ii], &hsd[ii], ptr_memory);
//       ptr_memory += (&hsd[ii])->memory_size;
//       d_cvt_vec2strvec(nb[ii], hlb[ii], &hsd[ii], 0);
//       d_cvt_vec2strvec(nb[ii], hub[ii], &hsd[ii], nb[ii]);
//       d_cvt_vec2strvec(ng[ii], hlg[ii], &hsd[ii], 2*nb[ii]);
//       d_cvt_vec2strvec(ng[ii], hug[ii], &hsd[ii], 2*nb[ii] + ng[ii]);
//
//       d_create_strvec(nu[ii]+nx[ii], &hsux[ii], ptr_memory);
//       ptr_memory += (&hsux[ii])->memory_size;
//
//       d_create_strvec(nx[ii+1], &hspi[ii], ptr_memory);
//       ptr_memory += (&hspi[ii])->memory_size;
//
//       d_create_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
//       ptr_memory += (&hslam[ii])->memory_size;
//
//       d_create_strvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
//       ptr_memory += (&hst[ii])->memory_size;
//     }
//
//     ii = N;
//
//     d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], ptr_memory);
//     ptr_memory += (&hsRSQrq[ii])->memory_size;
//     d_cvt_mat2strmat(nu[ii], nu[ii], hR[ii], nu[ii], &hsRSQrq[ii], 0, 0);
//     d_cvt_tran_mat2strmat(nu[ii], nx[ii], hS[ii], nu[ii], &hsRSQrq[ii], nu[ii], 0);
//     d_cvt_mat2strmat(nx[ii], nx[ii], hQ[ii], nx[ii], &hsRSQrq[ii], nu[ii], nu[ii]);
//     d_cvt_tran_mat2strmat(nu[ii], 1, hr[ii], nu[ii], &hsRSQrq[ii], nu[ii]+nx[ii], 0);
//     d_cvt_tran_mat2strmat(nx[ii], 1, hq[ii], nx[ii], &hsRSQrq[ii], nu[ii]+nx[ii], nu[ii]);
//
//     d_create_strvec(2*nb[ii]+2*ng[ii], &hsd[ii], ptr_memory);
//     ptr_memory += (&hsd[ii])->memory_size;
//     d_cvt_vec2strvec(nb[ii], hlb[ii], &hsd[ii], 0);
//     d_cvt_vec2strvec(nb[ii], hub[ii], &hsd[ii], nb[ii]);
//     d_cvt_vec2strvec(ng[ii], hlg[ii], &hsd[ii], 2*nb[ii]);
//     d_cvt_vec2strvec(ng[ii], hug[ii], &hsd[ii], 2*nb[ii] + ng[ii]);
//
//     d_create_strvec(nu[ii]+nx[ii], &hsux[ii], ptr_memory);
//     ptr_memory += (&hsux[ii])->memory_size;
//
//     d_create_strvec(nx[ii], &hspi[ii], ptr_memory);  // TODO(Andrea): bug?
//     ptr_memory += (&hspi[ii])->memory_size;
//
//     d_create_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
//     ptr_memory += (&hslam[ii])->memory_size;
//
//     d_create_strvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
//     ptr_memory += (&hst[ii])->memory_size;
//
//     // TODO(Andrea): LAST STAGE?
//
//     // extract output struct members
//     double **hx = qp_out->x;
//     double **hu = qp_out->u;
//     // double **hpi = qp_out->pi; //TODO(Andrea): multiplers not returned at the moment
//     // double **hlam = qp_out->lam; //TODO(Andrea): multiplers not returned at the moment
//
//     // extract args struct members
//     double mu_tol = hpmpc_args->tol;
//     int k_max = hpmpc_args->max_iter;
//     double mu0 = hpmpc_args->mu0;
//     int warm_start = hpmpc_args->warm_start;
//     // int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem
//
//     //  other solver arguments
//     int kk = -1;  // actual number of iterations
//     // double inf_norm_res[4];  // inf norm of residuals
//     // for (ii = 0; ii < 4; ii++) inf_norm_res[ii] = 0.0;  // zero
//
//     // double mu;
//     // struct d_strvec hswork[2];
//
//     // max sizes
//     int ngM = 0;
//     for ( ii=0; ii <= N; ii++ ) {
//       ngM = ng[ii] > ngM ? ng[ii] : ngM;
//     }
//
//     int nzM  = 0;
//     for ( ii =  0; ii <= N; ii++ ) {
//       nzM = nu[ii]+nx[ii]+1 > nzM ? nu[ii]+nx[ii]+1 : nzM;
//     }
//
//     int nxgM = ng[N];
//     for ( ii=0; ii < N; ii++ ) {
//       nxgM = nx[ii+1]+ng[ii] > nxgM ? nx[ii+1]+ng[ii] : nxgM;
//     }
//
//     // IPM constants
//     int hpmpc_status;
//     double alpha_min = 1e-8;
//     double *stat; d_zeros(&stat, k_max, 5);
//     // int compute_res = 0;
//     int compute_mult = 1;
//
//     // void *work_memory;
//     // v_zeros_align(&work_memory,
//       // d_ip2_res_mpc_hard_tv_work_space_size_bytes_libstr(N, nx, nu, nb, ng));
//
//     hpmpc_status = d_ip2_res_mpc_hard_libstr(&kk, k_max, mu0, mu_tol, alpha_min,
//       warm_start, stat, N, nx, nu, nb, hsidxb, ng, hsBAbt, hsRSQrq, hsDCt,
//       hsd, hsux, compute_mult, hspi, hslam, hst, ptr_memory);
//
//     double **temp_u;
//     // copy result to qp_out
//     for ( ii = 0; ii <= N; ii++ ) {
//       hu[ii] = hsux[ii].pa;
//       temp_u = &hsux[ii].pa;
//       hx[ii] = &temp_u[0][nu[ii]];
//     }
//
//     if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;
//
//     if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;
//
//     hpmpc_args->out_iter = kk;
//     // return
//     return acados_status;
// }

// int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
  // ocp_qp_hpmpc_args *hpmpc_args, void *workspace_) {
int ocp_qp_hpmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
        void *args_, void *mem_, void *workspace_) {
//
    ocp_qp_hpmpc_args *hpmpc_args = (ocp_qp_hpmpc_args*) args_;

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

    // extract output struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;  // TODO(Andrea): not returning multiplers atm
    double **hlam = qp_out->lam;
    double **ht = qp_out->t;

    int hpmpc_status = -1;




    int_t M = hpmpc_args->M;

    if (M < N) {  // XXX andrea partial tightening stuff
//
    // Process arguments TODO(Andrea): ask dimitris what this is for
    // args_args->dummy = 1.0;
    // workspace_++;
    // workspace_ = 0;
    // workspace_ = 0;
    // workspace_++;


    char *ptr_memory = (char *) mem_;


    // extract args struct members
    double mu_tol = hpmpc_args->tol;
    int k_max = hpmpc_args->max_iter;
    double mu0 = hpmpc_args->mu0;
    int warm_start = hpmpc_args->warm_start;

    //  other solver arguments
    int kk = -1;  // actual number of iterations

    // IPM constants
    // int kk_avg;
    double alpha_min = 1e-8;
    double *stat = (double*)ptr_memory;
    ptr_memory+=sizeof(double)*k_max*6;
    int compute_mult = 1;


    struct d_strmat *hsmatdummy = NULL;
    struct d_strvec *hsvecdummy = NULL;

    struct d_strmat hsBAbt[N+1];
    struct d_strvec hsb[N+1];
    struct d_strmat hsRSQrq[N+1];
    struct d_strvec hsQx[N+1];
    struct d_strvec hsqx[N+1];
    struct d_strvec hstinv[N+1];
    struct d_strvec hsrq[N+1];
    struct d_strmat hsDCt[N+1];
    struct d_strvec hsd[N+1];
    struct d_strvec hsux[N+1];
    struct d_strvec hsdux[N+1];
    struct d_strvec hspi[N+1];
    struct d_strvec hslam[N+1];
    struct d_strvec hst[N+1];
    struct d_strvec hsPb[N+1];
    struct d_strmat hsL[N+1];
//    struct d_strmat hsLxt[N+1];
    struct d_strmat hsric_work_mat[2];

    struct d_strvec hsdlam[N+1];  // to be checked
    struct d_strvec hsdt[N+1];
    struct d_strvec hslamt[N+1];  // to be checked

    for ( ii = 0; ii < N; ii++ ) {
      d_create_strmat(nu[ii]+nx[ii]+1, nx[ii+1], &hsBAbt[ii], ptr_memory);
      ptr_memory += (&hsBAbt[ii])->memory_size;
      d_cvt_tran_mat2strmat(nx[ii+1], nu[ii], hB[ii], nx[ii+1], &hsBAbt[ii], 0, 0);
      d_cvt_tran_mat2strmat(nx[ii+1], nx[ii], hA[ii], nx[ii+1], &hsBAbt[ii], nu[ii], 0);
      d_cvt_tran_mat2strmat(nx[ii+1], 1, hb[ii], nx[ii+1], &hsBAbt[ii], nu[ii]+nx[ii], 0);

      d_create_strvec(nx[ii+1], &hsb[ii], ptr_memory);
      ptr_memory += (&hsb[ii])->memory_size;
      d_cvt_vec2strvec(nx[ii+1], hb[ii], &hsb[ii], 0);

      d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], ptr_memory);
      ptr_memory += (&hsRSQrq[ii])->memory_size;
      d_cvt_mat2strmat(nu[ii], nu[ii], hR[ii], nu[ii], &hsRSQrq[ii], 0, 0);
      d_cvt_tran_mat2strmat(nu[ii], nx[ii], hS[ii], nu[ii], &hsRSQrq[ii], nu[ii], 0);
      d_cvt_mat2strmat(nx[ii], nx[ii], hQ[ii], nx[ii], &hsRSQrq[ii], nu[ii], nu[ii]);
      d_cvt_tran_mat2strmat(nu[ii], 1, hr[ii], nu[ii], &hsRSQrq[ii], nu[ii]+nx[ii], 0);
      d_cvt_tran_mat2strmat(nx[ii], 1, hq[ii], nx[ii], &hsRSQrq[ii], nu[ii]+nx[ii], nu[ii]);

      d_create_strvec(nu[ii]+nx[ii], &hsrq[ii], ptr_memory);
      ptr_memory += (&hsrq[ii])->memory_size;
      d_cvt_vec2strvec(nu[ii], hr[ii], &hsrq[ii], 0);
      d_cvt_vec2strvec(nx[ii], hq[ii], &hsrq[ii], nu[ii]);

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

      // initialize hsdux to primal input later usx will be subtracted
      d_create_strvec(nu[ii]+nx[ii], &hsdux[ii], ptr_memory);
      d_cvt_vec2strvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsdux[ii], 0);
      ptr_memory += (&hsdux[ii])->memory_size;
      d_create_strvec(nu[ii]+nx[ii], &hsux[ii], ptr_memory);
      d_cvt_vec2strvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsux[ii], 0);
      ptr_memory += (&hsux[ii])->memory_size;

      d_create_strvec(nx[ii+1], &hspi[ii], ptr_memory);
      ptr_memory += (&hspi[ii])->memory_size;

      d_create_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
      // copy multipliers from hpmpc_args
      d_cvt_vec2strvec(2*nb[ii]+2*ng[ii], hpmpc_args->lam0[ii], &hslam[ii], 0);
      ptr_memory += (&hslam[ii])->memory_size;

      d_create_strvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
      // copy slacks from hpmpc_args
      d_cvt_vec2strvec(2*nb[ii]+2*ng[ii], hpmpc_args->t0[ii], &hst[ii], 0);
      ptr_memory += (&hst[ii])->memory_size;

      d_create_strvec(nx[ii+1], &hsPb[ii+1], ptr_memory);
      ptr_memory += (&hsPb[ii+1])->memory_size;
      d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsL[ii], ptr_memory);
      ptr_memory += (&hsL[ii])->memory_size;
//      d_create_strmat(nx[ii], nx[ii], &hsLxt[ii], ptr_memory);
//      ptr_memory += (&hsLxt[ii])->memory_size;

      d_create_strvec(2*nb[ii]+2*ng[ii], &hstinv[ii], ptr_memory);
      ptr_memory += (&hstinv[ii])->memory_size;
      d_create_strvec(nb[ii]+ng[ii], &hsQx[ii], ptr_memory);
      ptr_memory += (&hsQx[ii])->memory_size;
      d_create_strvec(nb[ii]+ng[ii], &hsqx[ii], ptr_memory);
      ptr_memory += (&hsqx[ii])->memory_size;

      d_create_strvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], ptr_memory);
      ptr_memory += (&hsdlam[ii])->memory_size;
      d_create_strvec(2*nb[ii]+2*ng[ii], &hsdt[ii], ptr_memory);
      ptr_memory += (&hsdt[ii])->memory_size;
      d_create_strvec(2*nb[ii]+2*ng[ii], &hslamt[ii], ptr_memory);
      ptr_memory += (&hslamt[ii])->memory_size;
    }

    ii = N;
    d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsRSQrq[ii], ptr_memory);
    ptr_memory += (&hsRSQrq[ii])->memory_size;
    d_cvt_mat2strmat(nx[ii], nx[ii], hQ[ii], nx[ii], &hsRSQrq[ii], nu[ii], nu[ii]);
    d_cvt_tran_mat2strmat(nx[ii], 1, hq[ii], nx[ii], &hsRSQrq[ii], nu[ii]+nx[ii], nu[ii]);

    d_create_strvec(nu[ii]+nx[ii], &hsrq[ii], ptr_memory);
    ptr_memory += (&hsrq[ii])->memory_size;
    d_cvt_vec2strvec(nx[ii], hq[ii], &hsrq[ii], nu[ii]);

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

    // initialize hsdux to primal input later usx will be subtracted
    d_create_strvec(nu[ii]+nx[ii], &hsdux[ii], ptr_memory);
    d_cvt_vec2strvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsdux[ii], 0);
    ptr_memory += (&hsdux[ii])->memory_size;
    d_create_strvec(nu[ii]+nx[ii], &hsux[ii], ptr_memory);
    d_cvt_vec2strvec(nu[ii]+nx[ii], hpmpc_args->ux0[ii], &hsux[ii], 0);
    ptr_memory += (&hsux[ii])->memory_size;

    d_create_strvec(nx[ii], &hspi[ii], ptr_memory);  // Andrea: bug?
    ptr_memory += (&hspi[ii])->memory_size;

    d_create_strvec(2*nb[ii]+2*ng[ii], &hslam[ii], ptr_memory);
    d_cvt_vec2strvec(2*nb[ii]+2*ng[ii], hpmpc_args->lam0[ii], &hslam[ii], 0);
    ptr_memory += (&hslam[ii])->memory_size;

    d_create_strvec(2*nb[ii]+2*ng[ii], &hst[ii], ptr_memory);
    d_cvt_vec2strvec(2*nb[ii]+2*ng[ii], hpmpc_args->t0[ii], &hst[ii], 0);
    ptr_memory += (&hst[ii])->memory_size;

    d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], &hsL[ii], ptr_memory);
    ptr_memory += (&hsL[ii])->memory_size;
//    d_create_strmat(nx[ii], nx[ii], &hsLxt[ii], ptr_memory);
//    ptr_memory += (&hsLxt[ii])->memory_size;

    d_create_strvec(2*nb[ii]+2*ng[ii], &hslamt[ii], ptr_memory);
    ptr_memory += (&hslamt[ii])->memory_size;

    d_create_strvec(2*nb[ii]+2*ng[ii], &hstinv[ii], ptr_memory);
    ptr_memory += (&hstinv[ii])->memory_size;
    d_create_strvec(nb[ii]+ng[ii], &hsQx[ii], ptr_memory);
    ptr_memory += (&hsQx[ii])->memory_size;
    d_create_strvec(nb[ii]+ng[ii], &hsqx[ii], ptr_memory);
    ptr_memory += (&hsqx[ii])->memory_size;
    d_create_strvec(2*nb[ii]+2*ng[ii], &hsdlam[ii], ptr_memory);
    ptr_memory += (&hsdlam[ii])->memory_size;
    d_create_strvec(2*nb[ii]+2*ng[ii], &hsdt[ii], ptr_memory);
    ptr_memory += (&hsdt[ii])->memory_size;





    real_t sigma_mu = hpmpc_args->sigma_mu;

    int nuM;
    int nbM;

    struct d_strmat sLxM;
    struct d_strmat sPpM;

    d_create_strmat(nx[M]+1, nx[M], &sLxM, ptr_memory);
    ptr_memory += (&sLxM)->memory_size;
    d_create_strmat(nx[M]+1, nx[M], &sPpM, ptr_memory);
    ptr_memory += (&sPpM)->memory_size;

    struct d_strmat hstmpmat0;

    // riccati work space
    void *work_ric;

    work_ric = ptr_memory;
    ptr_memory+=d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);

    // update cost function matrices and vectors (box constraints)
    d_update_hessian_gradient_mpc_hard_libstr(N-M, &nx[M], &nu[M], &nb[M], &ng[M], \
      &hsd[M], sigma_mu, &hst[M], &hstinv[M], &hslam[M], &hslamt[M], &hsdlam[M], \
      &hsQx[M], &hsqx[M]);

    // backward riccati factorization and solution at the end
    d_back_ric_rec_sv_back_libstr(N-M, &nx[M], &nu[M], &nb[M], &hsidxb[M], &ng[M], \
      0, &hsBAbt[M], hsvecdummy, 1, &hsRSQrq[M], &hsrq[M], &hsDCt[M], &hsQx[M], \
      &hsqx[M], &hsux[M], 1, &hspi[M],  1, &hsPb[M], &hsL[M], work_ric);

    // extract chol factor of [P p; p' *]
    // TODO(Andrea): have m and n !!!!!
    dtrcp_l_libstr(nx[M], 1.0, &hsL[M], nu[M], nu[M], &sLxM, 0, 0);
    dgecp_libstr(1, nx[M], 1.0, &hsL[M], nu[M]+nx[M], nu[M], &sLxM, nx[M], 0);

    // d_print_strmat(nx[M]+1, nx[M], &sLxM, 0, 0);

    // recover [P p; p' *]
    dsyrk_ln_mn_libstr(nx[M]+1, nx[M], nx[M], 1.0, &sLxM, 0, 0, &sLxM, 0, 0, 0.0,
      &sPpM, 0, 0, &sPpM, 0, 0);

    // backup stage M
    nuM = nu[M];
    nbM = nb[M];
    hstmpmat0 = hsRSQrq[M];

    // update new terminal cost
    nu[M] = 0;
    nb[M] = 0;
    hsRSQrq[M] = sPpM;
    hsux[M].pa += nuM;

    // IPM at the beginning
    hpmpc_status = d_ip2_res_mpc_hard_libstr(&kk, k_max, mu0, mu_tol, alpha_min,
      warm_start, stat, M, nx, nu, nb, hsidxb, ng, hsBAbt, hsRSQrq, hsDCt,
      hsd, hsux, compute_mult, hspi, hslam, hst, ptr_memory);  // recover original stage M

    nu[M] = nuM;
    nb[M] = nbM;
    hsRSQrq[M] = hstmpmat0;
    hsux[M].pa -= nuM;

    // forward riccati solution at the end
    d_back_ric_rec_sv_forw_libstr(N-M, &nx[M], &nu[M], &nb[M], &hsidxb[M], &ng[M],
      0, &hsBAbt[M], hsvecdummy, 1, &hsRSQrq[M], &hsrq[M], hsmatdummy,
      &hsQx[M], &hsqx[M], &hsux[M], 1, &hspi[M], 1, &hsPb[M], &hsL[M],
      hsric_work_mat);

    // compute alpha, dlam and dt
    real_t alpha = 1.0;
    // compute primal step hsdux for stages M to N
    real_t *temp_p1, *temp_p2;
    for (int_t i = M; i <= N; i++) {
      // hsdux is initialized to be equal to hpmpc_args.ux0
      temp_p1 = hsdux[i].pa;
      temp_p2 = hsux[i].pa;  // hsux[i].pa;
      for (int_t j = 0; j < nx[i]+nu[i]; j++) temp_p1[j]= - temp_p1[j] + temp_p2[j];
    }

    d_compute_alpha_mpc_hard_libstr(N-M, &nx[M], &nu[M], &nb[M], &hsidxb[M],
      &ng[M], &alpha, &hst[M], &hsdt[M], &hslam[M], &hsdlam[M], &hslamt[M],
      &hsdux[M], &hsDCt[M], &hsd[M]);

    // overwrite alpha (taking full steps and performing line-search in out_iter
    // level)

    // alpha = 1.0;

    // update stages M to N
    double mu_scal = 0.0;
    d_update_var_mpc_hard_libstr(N-M, &nx[M], &nu[M], &nb[M], &ng[M],
      &sigma_mu, mu_scal, alpha, &hsux[M], &hsdux[M], &hst[M], &hsdt[M], &hslam[M],
      &hsdlam[M], &hspi[M], &hspi[M]);

    // !!!! TODO(Andrea): equality multipliers are not being updated! Need to
    // define and compute hsdpi (see function prototype).

//    double **temp_u;
    // copy result to qp_out
    for ( ii = 0; ii < N; ii++ ) {
        d_cvt_strvec2vec(nu[ii], &hsux[ii], 0, hu[ii]);
        d_cvt_strvec2vec(nx[ii], &hsux[ii], nu[ii], hx[ii]);
        d_cvt_strvec2vec(nx[ii], &hspi[ii], 0, hpi[ii]);
        d_cvt_strvec2vec(2*nb[ii]+2*ng[ii], &hslam[ii], 0, hlam[ii]);
        d_cvt_strvec2vec(2*nb[ii]+2*ng[ii], &hst[ii], 0, ht[ii]);
//      hu[ii] = hsux[ii].pa;
//      hlam[ii] = hslam[ii].pa;
//      ht[ii] = hst[ii].pa;
//      temp_u = &hsux[ii].pa;
//      hx[ii] = &temp_u[0][nu[ii]];
    }

    ii = N;
    d_cvt_strvec2vec(nx[ii], &hsux[ii], nu[ii], hx[ii]);
    d_cvt_strvec2vec(2*nb[ii]+2*ng[ii], &hslam[ii], 0, hlam[ii]);
    d_cvt_strvec2vec(2*nb[ii]+2*ng[ii], &hst[ii], 0, ht[ii]);
//    hlam[ii] = hslam[ii].pa;
//    ht[ii] = hst[ii].pa;
//    temp_u = &hsux[ii].pa;
//    hx[ii] = &temp_u[0][nu[ii]];

    hpmpc_args->out_iter = kk;
//
    } else {  // XXX giaf fortran order interface with partial condensing
//
     // extract args struct members
     double mu_tol = hpmpc_args->tol;
     int k_max = hpmpc_args->max_iter;
     double mu0 = hpmpc_args->mu0;
     int warm_start = hpmpc_args->warm_start;
     int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem
     int out_iter = -1;  // number of performed iterations
     double *inf_norm_res = hpmpc_args->inf_norm_res;

     // memory for stat
     size_t addr = ( ( (size_t) workspace_) + 7) / 8 * 8;  // align to 8-byte boundaries
     double *ptr_double = (double *) addr;
     double *stat = ptr_double;
     ptr_double += 5*k_max;
     for (ii = 0; ii < 5*k_max; ii++) stat[ii] = 0.0;  // zero

     // memory for idxb
     int *hidxb[N+1];
//   size_t addr = (( (size_t) workspace_ ) + 3 ) / 4 * 4;  // align to 4-byte boundaries
     int *ptr_int = (int *) ptr_double;
     for (ii = 0; ii <= N; ii++) {
         hidxb[ii] = ptr_int;
         ptr_int += nb[ii];
     }
     void *ipm_workspace = (void *) ptr_int;

     //  swap x and u in bounds (by updating their indeces)
     for (ii = 0; ii <= N; ii++) {
         jj = 0;
         for (; jj < nb[ii]; jj++) {
             if (hsidxb[ii][jj] < nx[ii]) {  // state
                 hidxb[ii][jj] = hsidxb[ii][jj]+nu[ii];
             } else {  // input
                 hidxb[ii][jj] = hsidxb[ii][jj]-nx[ii];
             }
         }
     }

     hpmpc_status = fortran_order_d_ip_ocp_hard_tv(&out_iter, k_max, mu0,
        mu_tol, N, nx, nu, nb, hidxb, ng, N2, warm_start, hA, hB, hb, hQ, hS,
        hR, hq, hr, hlb, hub, hC, hD, hlg,
         hug, hx, hu, hpi, hlam, inf_norm_res, ipm_workspace, stat);

     hpmpc_args->out_iter = out_iter;  // number of performed iterations
//
    }

    if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;
    if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;

    // return
    return acados_status;
}





int_t ocp_qp_hpmpc_calculate_workspace_size(ocp_qp_in *qp_in, void *args_) {
//
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args*) args_;

    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int **hidxb = (int **) qp_in->idxb;
    int *ng = (int *) qp_in->nc;
    int_t N2 = args->N2;
    int_t M = args->M;

    int_t ws_size;  // TODO(Andrea): dummy expression. Need to
    // ws_size = 0*in->N*args->N;  // TODO(Andrea): dummy expression. Need to
    // decide what hpmpc's workpsace is.

    if (M < N) {  // XXX andrea partial tightening stuff
        ws_size = 0;
    } else {  // XXX giaf fortran interface
        int ii;
        int_t max_ip_iter = args->max_iter;
        ws_size = 8 + 5*max_ip_iter*sizeof(double);
        for (ii = 0; ii <= N; ii++) {
            ws_size += nb[ii]*sizeof(int);
        }
        ws_size += hpmpc_d_ip_ocp_hard_tv_work_space_size_bytes(N, nx, \
            nu, nb, hidxb, ng, N2);
    }
    return ws_size;
}





int_t ocp_qp_hpmpc_create_memory(ocp_qp_in *in, void *args_, void **mem_) {
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args*) args_;
    // ocp_qp_hpmpc_memory mem = (ocp_qp_hpmpc_memory *) mem_;

    int_t N = (int_t)in->N;
    int_t *nx = (int_t*)in->nx;
    int_t *nu = (int_t*)in->nu;
    int_t *nb = (int_t*)in->nb;
//    int **hidxb = (int **) in->idxb;
    int_t *ng = (int_t*)in->nc;

    int_t max_ip_iter = args->max_iter;

    int_t M = args->M;

    int_t mem_size;

    if (M < N) {  // XXX andrea partial tightened stuff
//
    mem_size = d_ip2_res_mpc_hard_work_space_size_bytes_libstr(N,
      nx, nu, nb, ng);

    int_t ii;
    // Adding memory for data
    for (ii=0; ii < N; ii++) {
      mem_size+= d_size_strmat(nu[ii]+nx[ii]+1, nx[ii+1]);  // BAbt
      mem_size+= d_size_strvec(nx[ii+1]);  // b
      mem_size+= d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]);  // RSQrq
      mem_size+= d_size_strvec(nu[ii]+nx[ii]);  // rq
      mem_size+= d_size_strmat(nu[ii]+nx[ii]+1, ng[ii]);  // DCt
      mem_size+= d_size_strvec(2*nb[ii]+2*ng[ii]);  // d
      mem_size+= d_size_strvec(nu[ii]+nx[ii]);
      mem_size+= d_size_strvec(nx[ii+1]);
      mem_size+= d_size_strvec(nx[ii+1]);
      mem_size+= d_size_strvec(2*nb[ii]+2*ng[ii]);
      mem_size+= d_size_strvec(2*nb[ii]+2*ng[ii]);
      mem_size+= d_size_strvec(2*nb[ii]+2*ng[ii]);
    }

    mem_size+= d_size_strmat(nu[N]+nx[N]+1, nu[N]+nx[N]);  // RSQrq
    mem_size+= d_size_strvec(nu[N]+nx[N]);  // q
    mem_size+= d_size_strmat(nu[N]+nx[N]+1, ng[N]);  // DCt
    mem_size+= d_size_strvec(2*nb[N]+2*ng[N]);  // d
    mem_size+= d_size_strvec(nu[N]+nx[N]);
    mem_size+= d_size_strvec(nu[N]+nx[N]);
    mem_size+= d_size_strvec(2*nb[N]+2*ng[N]);
    mem_size+= d_size_strvec(2*nb[N]+2*ng[N]);
    mem_size+= d_size_strvec(2*nb[N]+2*ng[N]);


    // Adding memory for extra variables in the Riccati recursion
    mem_size+=d_size_strvec(nx[ii+1]);
    mem_size+=d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]);
    mem_size+=d_size_strmat(nx[ii], nx[ii]);

    for ( int ii=0; ii < N; ii++ ) {
      mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);
      mem_size+=d_size_strvec(nb[ii]+ng[ii]);
      mem_size+=d_size_strvec(nb[ii]+ng[ii]);

      mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);
      mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);
    }

    mem_size+=d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]);
    mem_size+=d_size_strmat(nx[ii], nx[ii]);

    ii = N;
    mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);
    mem_size+=d_size_strvec(nb[ii]+ng[ii]);
    mem_size+=d_size_strvec(nb[ii]+ng[ii]);
    mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);
    mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);

    mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);
    mem_size+=d_size_strvec(nb[ii]+ng[ii]);
    mem_size+=d_size_strvec(nb[ii]+ng[ii]);
    mem_size+=d_size_strvec(2*nb[ii]+2*ng[ii]);

    // mem_size+=d_size_strmat(nx[M]+1, nx[M]);
    // mem_size+=d_size_strmat(nx[M]+1, nx[M]);


    mem_size+=d_back_ric_rec_work_space_size_bytes_libstr(N, nx, nu, nb, ng);
    // add memory for riccati work space

    mem_size+=sizeof(double)*max_ip_iter*5;
    mem_size+=sizeof(double)*1000000;
    // add memory for stats

    } else {  // XXX giaf fortran interface
        mem_size = 0;
    }

    v_zeros_align(mem_, mem_size);

    return 1;
  }

void ocp_qp_hpmpc_free_memory(void *mem_) {
    ocp_qp_hpmpc_memory *mem = (ocp_qp_hpmpc_memory *) mem_;
    free(mem);
}

int_t ocp_qp_hpmpc_create_arguments(void *args_, int_t opts_) {
    ocp_qp_hpmpc_args *args = (ocp_qp_hpmpc_args*) args_;
    hpmpc_options_t opts = (hpmpc_options_t) opts_;

    if (opts == HPMPC_DEFAULT_ARGUMENTS) {
    args->tol = 1e-8;
    args->max_iter = 20;
    args->mu0 = 0.1;
    args->warm_start = 0;
    args->M = args->N;
    } else {
      printf("Invalid hpmpc options.");
      return -1;
    }
  return 0;
  }


// TODO(Andrea): need to merge hpmpc in order to use this...
// int ocp_qp_hpnmpc(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_hpmpc_args *hpmpc_args,
//     void *workspace_) {
//
//     // initialize return code
//     int acados_status = ACADOS_SUCCESS;
//
//     // loop index
//     int ii, jj;
//
//     // extract input struct members
//     int N = qp_in->N;
//     int *nx = (int *) qp_in->nx;
//     int *nu = (int *) qp_in->nu;
//     int *nb = (int *) qp_in->nb;
//     int **hidxb_swp = (int **) qp_in->idxb;
//     int *ng = (int *) qp_in->nc;
//     double **hA = (double **) qp_in->A;
//     double **hB = (double **) qp_in->B;
//     double **hb = (double **) qp_in->b;
//     double **hQ = (double **) qp_in->Q;
//     double **hS = (double **) qp_in->S;
//     double **hR = (double **) qp_in->R;
//     double **hq = (double **) qp_in->q;
//     double **hr = (double **) qp_in->r;
//     double **hlb = (double **) qp_in->lb;
//     double **hub = (double **) qp_in->ub;
//     double **hC = (double **) qp_in->Cx;
//     double **hD = (double **) qp_in->Cu;
//     double **hlg = (double **) qp_in->lc;
//     double **hug = (double **) qp_in->uc;
//
//     // extract output struct members
//     double **hx = qp_out->x;
//     double **hu = qp_out->u;
//     double **hpi = qp_out->pi;
//     double **hlam = qp_out->lam;
//     double **ht = qp_out->t;
//
//     // extract args struct members
//     // double mu_tol = hpmpc_args->tol;
//     int k_max = hpmpc_args->max_iter;
//     double mu0 = hpmpc_args->mu0;
//     int warm_start = hpmpc_args->warm_start;
//     int N2 = hpmpc_args->N2;  // horizon length of the partially condensed problem
//
//     double **ux0 = hpmpc_args->ux0;
//     double **pi0 = hpmpc_args->pi0;
//     double **lam0 = hpmpc_args->lam0;
//     double **t0 = hpmpc_args->t0;
//
//     //  other solver arguments
//     int kk = -1;  // actual number of iterations
//     // double inf_norm_res[4];  // inf norm of residuals
//     // for (ii = 0; ii < 4; ii++) inf_norm_res[ii] = 0.0;  // zero
//
//     // memory for stat
//     size_t addr = ( ( (size_t) workspace) + 7) / 8 * 8;  // align to 8-byte boundaries
//     double *ptr_double = (double *) addr;
//     double *stat = ptr_double;
//     ptr_double += 5*k_max;
//     for (ii = 0; ii < 5*k_max; ii++) stat[ii] = 0.0;  // zero
//
//     // memory for idxb
//     int *hidxb[N+1];
// //  size_t addr = (( (size_t) workspace ) + 3 ) / 4 * 4;  // align to 4-byte boundaries
//     int *ptr_int = (int *) ptr_double;
//     for (ii = 0; ii <= N; ii++) {
//         hidxb[ii] = ptr_int;
//         ptr_int += nb[ii];
//     }
//     workspace = (void *) ptr_int;
//
//     //  swap x and u in bounds (by updating their indeces)
//     for (ii = 0; ii <= N; ii++) {
//         jj = 0;
//         for (; jj < nb[ii]; jj++) {
//             if (hidxb_swp[ii][jj] < nx[ii]) {  // state
//                 hidxb[ii][jj] = hidxb_swp[ii][jj]+nu[ii];
//             } else {  // input
//                 hidxb[ii][jj] = hidxb_swp[ii][jj]-nx[ii];
//             }
//         }
//     }
// //
//     int hpmpc_status;
//     // hpmpc_status = fortran_order_d_ip_ocp_hard_tv_single_newton_step(&kk,
//     // k_max, mu0, mu_tol, N, nx, nu, nb,
//     // hidxb, ng, N2, warm_start, hA, hB, hb, hQ, hS, hR, hq, hr, hlb, hub, hC, hD, hlg, hug,
//     // hx, hu, hpi, hlam, ht, inf_norm_res, workspace, stat, ux0, pi0, lam0, t0);
//
//     if (hpmpc_status == 1) acados_status = ACADOS_MAXITER;
//
//     if (hpmpc_status == 2) acados_status = ACADOS_MINSTEP;
//
//     // return
//     return acados_status;
// }

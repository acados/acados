/**************************************************************************************************
* acados/external/hpmpc/include/mpc_solvers.h                                                                                                *
* This file is part of HPMPC.                                                                     *
*                                                                                                 *
* HPMPC -- Library for High-Performance implementation of solvers for MPC.                        *
* Copyright (C) 2014-2015 by Technical University of Denmark. All rights reserved.                *
*                                                                                                 *
* HPMPC is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPMPC is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPMPC; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, giaf (at) dtu.dk                                                       *
*                                                                                                 *
**************************************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

// 2017.03.13 Dang add target.h
#include "target.h"


// IPM without residuals computation
int d_ip2_mpc_hard_tv_work_space_size_bytes(int N, int *nx, int *nu, int *nb, int *ng);
int d_ip2_mpc_hard_tv(int *kk, int k_max, double mu0, double mu_tol, double alpha_min, int warm_start, double *stat, int N, int *nx, int *nu, int *nb, int **idxb, int *ng, double **pBAbt, double **pQ, double **pDCt, double **d, double **ux, int compute_mult, double **pi, double **lam, double **t, double *double_work_memory);
void d_kkt_solve_new_rhs_mpc_hard_tv(int N, int *nx, int *nu_N, int *nb, int **idxb, int *ng, double **pBAbt, double **r_A, double **pQ, double **r_H, double **pDCt, double **r_C, double **ux, int compute_mult, double **pi, double **lam, double **t, double *double_work_memory);
void d_res_mpc_hard_tv(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, double **hpBAbt, double **hb, double **hpQ, double **hq, double **hux, double **hpDCt, double **hd, double **hpi, double **hlam, double **ht, double **hrq, double **hrb, double **hrd, double *mu);



// IPM with residuals computation
int d_ip2_res_mpc_hard_tv_work_space_size_bytes(int N, int *nx, int *nu, int *nb, int *ng);
int d_ip2_res_mpc_hard_tv(int *kk, int k_max, double mu0, double mu_tol, double alpha_min, int warm_start, double *stat, int N, int *nx, int *nu_N, int *nb, int **idxb, int *ng, double **pBAbt, double **pQ, double **pDCt, double **d, double **ux, int compute_mult, double **pi, double **lam, double **t, double *double_work_memory);
void d_kkt_solve_new_rhs_res_mpc_hard_tv(int N, int *nx, int *nu_N, int *nb, int **idxb, int *ng, double **pBAbt, double **b, double **pQ, double **q, double **pDCt, double **d, double **ux, int compute_mult, double **pi, double **lam, double **t, double *double_work_memory);
//void d_res_res_mpc_hard_tv(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, double **hpBAbt, double **hb, double **hpQ, double **hq, double **hux, double **hpDCt, double **hd, double **hpi, double **hlam, double **ht, double *work, double **hrq, double **hrb, double **hrd, double **hrm, double *mu);



#ifdef BLASFEO
int d_ip2_res_mpc_hard_work_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int *ng);
int d_ip2_res_mpc_hard_libstr(int *kk, int k_max, double mu0, double mu_tol, double alpha_min, int warm_start, double *stat, int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsd, struct d_strvec *hsux, int compute_mult, struct d_strvec *hspi, struct d_strvec *hslam, struct d_strvec *hst, void *work_memory);
//int d_res_res_mpc_hard_work_space_size_bytes_libstr(int N, int *nx, int *nu, int *nb, int *ng);
//void d_res_res_mpc_hard_libstr(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct d_strmat *hsBAbt, struct d_strvec *hsb, struct d_strmat *hsQ, struct d_strvec *hsq, struct d_strvec *hsux, struct d_strmat *hsDCt, struct d_strvec *hsd, struct d_strvec *hspi, struct d_strvec *hslam, struct d_strvec *hst, struct d_strvec *hsrq, struct d_strvec *hsrb, struct d_strvec *hsrd, struct d_strvec *hsrm, double *mu, void *work);
#endif
#if defined(TREE_MPC)
#ifdef BLASFEO
int d_tree_ip2_res_mpc_hard_work_space_size_bytes_libstr(int Nn, struct node *tree, int *nx, int *nu, int *nb, int *ng);
int d_tree_ip2_res_mpc_hard_libstr(int *kk, int k_max, double mu0, double mu_tol, double alpha_min, int warm_start, double *stat, int Nn, struct node *tree, int *nx, int *nu, int *nb, int **idxb, int *ng, struct d_strmat *hsBAbt, struct d_strmat *hsRSQrq, struct d_strmat *hsDCt, struct d_strvec *hsd, struct d_strvec *hsux, int compute_mult, struct d_strvec *hspi, struct d_strvec *hslam, struct d_strvec *hst, void *work_memory);
#endif
#endif


// soft constraints // TODO old version
int d_ip2_mpc_soft_tv_work_space_size_bytes(int N, int *nx, int *nu, int *nb, int *ng, int *ns);
int d_ip2_mpc_soft_tv(int *kk, int k_max, double mu0, double mu_tol, double alpha_min, int warm_start, double *stat, int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int *ns, double **pBAbt, double **pQ, double **Z, double **z, double **pDCt, double **d, double **ux, int compute_mult, double **pi, double **lam, double **t, double *double_work_memory);
void d_res_mpc_soft_tv(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int *ns, double **hpBAbt, double **hpQ, double **hq, double **hZ, double **hz, double **hux, double **hpDCt, double **hd, double **hpi, double **hlam, double **ht, double **hrq, double **hrb, double **hrd, double **hrz, double *mu);



#ifdef __cplusplus
}
#endif

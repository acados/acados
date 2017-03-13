/**************************************************************************************************
* /acados/external/hpmpc/include/lqcp_aux.h                                                                                                *
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



//void d_update_diag_pmat(int kmax, double *pQ, int sda, double *d);
//void d_update_diag_pmat_sparse(int kmax, int *idx, double *pQ, int sda, double *d);
//void d_update_row_pmat(int kmax, double *pQ, double *r);
//void d_update_row_pmat_sparse(int kmax, int *idx, double *pQ, double *r);
void d_update_vector_sparse(int kmax, int *idx, double *q, double *v);
//void d_add_row_pmat(int kmax, double *pA, double *pC);
//void d_add_diag_pmat(int kmax, double *pA, int sda, double *d);
void dpotrf_diag_lib(int m, int n, double *C, int sdc, double *D, int sdd);



#ifdef __cplusplus
}
#endif

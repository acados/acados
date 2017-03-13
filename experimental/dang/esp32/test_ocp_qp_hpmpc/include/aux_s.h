/**************************************************************************************************
* acados/external/hpmpc/include/aux_s.h                                                                                                *
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



void s_zeros(float **pA, int row, int col);
void s_zeros_align(float **pA, int row, int col);
void s_eye(float **pA, int row);
void s_copy_mat(int row, int col, float *A, int lda, float *B, int ldb);
void s_copy_pmat(int row, int col, int bs, float *A, int sda, float *B, int sdb);
//void s_copy_pmat_lo(int row, int bs, float *A, int sda, float *B, int sdb);
//void s_transpose_pmat_lo(int row, int offset, float *A, int sda, float *B, int sdb);
//void s_align_pmat(int row, int col, int offset, int bs, float *A, int sda, float *B, int sdb);
void s_cvt_d2s_pmat(int row, int col, int dbs, double *A, int sda, int sbs, float *B, int sdb);
void s_cvt_mat2pmat(int row, int col, float *A, int lda, int offset, float *B, int sdb);
void cvt_d2s_mat2pmat(int row, int col, int offset, int bs_dummy, double *A, int lda, float *pA, int sda);
void s_cvt_tran_mat2pmat(int row, int col, int offset, int bs, float *A, int lda, float *B, int sdb);
void cvt_tran_d2s_mat2pmat(int row, int col, int offset, int bs_dummy, double *A, int lda, float *pA, int sda);
void s_cvt_pmat2mat(int row, int col, int offset, int bs, float *pA, int sda, float *A, int lda);
void s_print_mat(int row, int col, float *A, int lda);
void s_print_pmat(int row, int col, int bs, float *A, int sda);



#ifdef __cplusplus
}
#endif

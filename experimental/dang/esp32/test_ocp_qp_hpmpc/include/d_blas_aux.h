/**************************************************************************************************
* acados/external/hpmpc/include/d_blas_aux.h                                                                                                *
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


// level 1 BLAS
#if ! defined(BLASFEO)
void daxpy_lib(int kmax, double alpha, double *x, double *y);
#endif
void daxpy_bkp_lib(int kmax, double alpha, double *x, double *y, double *z);



// auxiliary routines
#if ! defined(BLASFEO)
void ddiareg_lib(int kmax, double reg, int offset, double *pD, int sdd);
void ddiain_lib(int kmax, double *x, int offset, double *pD, int sdd);
void ddiain_sqrt_lib(int kmax, double *x, int offset, double *pD, int sdd);
void ddiaex_lib(int kmax, int offset, double *pD, int sdd, double *x);
void ddiaad_lib(int kmax, double alpha, double *x, int offset, double *pD, int sdd);
void ddiain_libsp(int kmax, int *idx, double *x, double *pD, int sdd);
void ddiaad_libsp(int kmax, int *idx, double alpha, double *x, double *pD, int sdd);
void ddiaadin_libsp(int kmax, int *idx, double alpha, double *x, double *y, double *pD, int sdd);
void drowin_lib(int kmax, double *x, double *pD);
void drowex_lib(int kmax, double *pD, double *x);
void drowad_lib(int kmax, double alpha, double *x, double *pD);
void drowin_libsp(int kmax, int *idx, double *x, double *pD);
void drowad_libsp(int kmax, int *idx, double alpha, double *x, double *pD);
void drowsw_lib(int kmax, double *pA, double *pC);
void dcolin_lib(int kmax, double *x, int offset, double *pD, int sdd);
void dcolad_lib(int kmax, double alpha, double *x, int offset, double *pD, int sdd);
void dcolin_libsp(int kmax, int *idx, double *x, double *pD, int sdd);
void dcolad_libsp(int kmax, double alpha, int *idx, double *x, double *pD, int sdd);
void dvecad_libsp(int kmax, int *idx, double alpha, double *x, double *y);
#endif




// diagonal routines
#if ! defined(BLASFEO)
void dgemm_diag_right_lib(int m, int n, double *pA, int sda, double *dB, int alg, double *pC, int sdc, double *pD, int sdd);
void dgemm_diag_left_lib(int m, int n, double *dA, double *pB, int sdb, int alg, double *pC, int sdc, double *pD, int sdd);
#endif
void dsyrk_diag_left_right_lib(int m, double *Al, double *Ar, double *B, int sdb, int alg, double *C, int sdc, double *D, int sdd);
void dgemv_diag_lib(int m, double *dA, double *x, int alg, double *y, double *z);



#ifdef __cplusplus
}
#endif

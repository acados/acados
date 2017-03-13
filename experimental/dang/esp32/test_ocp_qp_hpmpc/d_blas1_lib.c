/**************************************************************************************************
* acados/external/blasfeo/blas/d_blas1_lib.c                                                                                                *
* This file is part of BLASFEO.                                                                   *
*                                                                                                 *
* BLASFEO -- BLAS For Embedded Optimization.                                                      *
* Copyright (C) 2016-2017 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
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
*                          gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#if defined(LA_BLAS)
#if defined(LA_BLAS_OPENBLAS)
#include <f77blas.h>
#elif defined(LA_BLAS_BLIS)
#elif defined(LA_BLAS_NETLIB)
#include "d_blas.h"
#elif defined(LA_BLAS_MKL)
#include <mkl_blas.h>
#endif
#endif

#include "blasfeo_common.h"
#include "blasfeo_d_kernel.h"



#if defined(LA_REFERENCE)



void daxpy_libstr(int m, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi)
	{
	int ii;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	for(ii=0; ii<m; ii++)
		y[ii] += alpha * x[ii];
	return;
	}



void daxpy_bkp_libstr(int m, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	int ii;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
	for(ii=0; ii<m; ii++)
		{
		z[ii] = y[ii];
		y[ii] += alpha * x[ii];
		}
	return;
	}



#elif defined(LA_BLAS)



void daxpy_libstr(int m, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi)
	{
	int i1 = 1;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
#if defined(LA_BLAS_MKL)
	daxpy(&m, &alpha, x, &i1, y, &i1);
#else
	daxpy_(&m, &alpha, x, &i1, y, &i1);
#endif
	return;
	}



void daxpy_bkp_libstr(int m, double alpha, struct d_strvec *sx, int xi, struct d_strvec *sy, int yi, struct d_strvec *sz, int zi)
	{
	int i1 = 1;
	double *x = sx->pa + xi;
	double *y = sy->pa + yi;
	double *z = sz->pa + zi;
#if defined(LA_BLAS_MKL)
	dcopy(&m, y, &i1, z, &i1);
	daxpy(&m, &alpha, x, &i1, y, &i1);
#else
	dcopy_(&m, y, &i1, z, &i1);
	daxpy_(&m, &alpha, x, &i1, y, &i1);
#endif
	return;
	}



#else

#error : wrong LA choice

#endif

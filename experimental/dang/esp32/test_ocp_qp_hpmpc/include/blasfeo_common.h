/**************************************************************************************************
* acados/external/blasfeo/include/blasfeo_common.h                                                                                                *
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


#ifdef __cplusplus
extern "C" {
#endif

// Dang added for make
#include "blasfeo_target.h"


#if defined(LA_HIGH_PERFORMANCE)

// matrix structure
struct d_strmat
	{
	int m; // rows
	int n; // cols
	int pm; // packed number or rows
	int cn; // packed number or cols
	double *pA; // pointer to a pm*pn array of doubles, the first is aligned to cache line size
	double *dA; // pointer to a min(m,n) (or max???) array of doubles
	int use_dA; // flag to tell if dA can be used
	int memory_size; // size of needed memory
	};

// vector structure
struct d_strvec
	{
	int m; // size
	int pm; // packed size
	double *pa; // pointer to a pm array of doubles, the first is aligned to cache line size
	int memory_size; // size of needed memory
	};

#elif defined(LA_BLAS) | defined(LA_REFERENCE)

// matrix structure
struct d_strmat
	{
	int m; // rows
	int n; // cols
	double *pA; // pointer to a m*n array of doubles
#if defined(LA_REFERENCE)
	double *dA; // pointer to a min(m,n) (or max???) array of doubles
	int use_dA; // flag to tell if dA can be used
#endif
	int memory_size; // size of needed memory
	};

// vector structure
struct d_strvec
	{
	int m; // size
	double *pa; // pointer to a m array of doubles, the first is aligned to cache line size
	int memory_size; // size of needed memory
	};

#else

#error : wrong LA choice

#endif



#ifdef __cplusplus
}
#endif

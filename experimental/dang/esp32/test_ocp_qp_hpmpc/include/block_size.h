/**************************************************************************************************
* /home/dang/acados/external/hpmpc/include/block_size.h                                                                                                *
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
// Manually define the target for compilation without HPMPC options
// 2017.03.13 Dang add target.h
#include "target.h"


#if defined( TARGET_X64_AVX2 )

#define D_MR 4
#define S_MR 8
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 2
#endif

#elif defined( TARGET_X64_AVX )

#define D_MR 4
#define S_MR 8
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 2
#endif

#elif defined( TARGET_X64_SSE3 )

#define D_MR 4
#define S_MR 4
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 4
#endif

#elif defined( TARGET_C99_4X4 )

#define D_MR 4
#define S_MR 4
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 4
#endif

#elif defined( TARGET_C99_4X4_PREFETCH )

#define D_MR 4
#define S_MR 4
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 4
#endif

#elif defined( TARGET_CORTEX_A57 )

#define D_MR 4
#define S_MR 4
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 4
#endif

#elif defined( TARGET_CORTEX_A15 )

#define D_MR 4
#define S_MR 4
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 4
#endif

#elif defined( TARGET_CORTEX_A9 )

#define D_MR 4
#define S_MR 4
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2 // 1
#define S_NCL 2
#endif

#elif defined( TARGET_CORTEX_A7 )

#define D_MR 4
#define S_MR 4
#if defined BLASFEO // XXX
#define D_NCL 4
#define S_NCL 4
#else
#define D_NCL 2
#define S_NCL 4
#endif


#else
#error "Unknown architecture"
#endif

/**************************************************************************************************
* acados/external/blasfeo/include/blasfeo_block_size.h                                                                                                *
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
// 2017.03.13: Dang added blasfeo_target.h
#include "blasfeo_target.h"

#if defined( TARGET_X64_INTEL_HASWELL )

#define D_BS 4
#define S_BS 4 //8
#define D_NC 4 // 2 // until the smaller kernel is 4x4
#define S_NC 4 //2

#elif defined( TARGET_X64_INTEL_SANDY_BRIDGE )

#define D_BS 4
#define S_BS 4 //8
#define D_NC 4 // 2 // until the smaller kernel is 4x4
#define S_NC 4 //2

#elif defined( TARGET_X64_INTEL_CORE )

#define D_BS 4
#define S_BS 4
#define D_NC 4 // 2 // until the smaller kernel is 4x4
#define S_NC 4 //2

#elif defined( TARGET_X64_AMD_BULLDOZER )

#define D_BS 4
#define S_BS 4
#define D_NC 4 // 2 // until the smaller kernel is 4x4
#define S_NC 4 //2

#elif defined( TARGET_GENERIC )

#define D_BS 4
#define S_BS 4
#define D_NC 4 // 2 // until the smaller kernel is 4x4
#define S_NC 4 //2

#else
#error "Unknown architecture"
#endif

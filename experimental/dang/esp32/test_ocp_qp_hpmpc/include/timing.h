/*    /home/dang/acados/acados/utils/timing.h
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

#ifndef ACADOS_UTILS_TIMING_H_
#define ACADOS_UTILS_TIMING_H_

#if defined(__APPLE__)
#include <mach/mach_time.h>
#else
#include <sys/stat.h>
#include <sys/time.h>
#endif

#include "types.h"

#if defined(__APPLE__)
typedef struct acado_timer_ {
    uint64_t tic;
    uint64_t toc;
    mach_timebase_info_data_t tinfo;
} acado_timer;
#else
typedef struct acado_timer_ {
    struct timeval tic;
    struct timeval toc;
} acado_timer;
#endif

void acado_tic(acado_timer* t);
real_t acado_toc(acado_timer* t);

#endif  // ACADOS_UTILS_TIMING_H_

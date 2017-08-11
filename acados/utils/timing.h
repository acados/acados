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

#ifndef ACADOS_UTILS_TIMING_H_
#define ACADOS_UTILS_TIMING_H_

#include "acados/utils/types.h"

#if !(defined __DSPACE__)
#if (defined _WIN32 || defined _WIN64) && !(defined __MINGW32__ || defined __MINGW64__)

/* Use Windows QueryPerformanceCounter for timing. */
#include <Windows.h>

/** A structure for keeping internal timer data. */
typedef struct acados_timer_ {
    LARGE_INTEGER tic;
    LARGE_INTEGER toc;
    LARGE_INTEGER freq;
} acados_timer;


#elif(defined __APPLE__)

#include <mach/mach_time.h>

/** A structure for keeping internal timer data. */
typedef struct acados_timer_ {
    uint64_t tic;
    uint64_t toc;
    mach_timebase_info_data_t tinfo;
} acados_timer;

#else

/* Use POSIX clock_gettime() for timing on non-Windows machines. */
#include <time.h>

#if __STDC_VERSION__ >= 199901L
/* C99 mode of operation. */

#include <sys/stat.h>
#include <sys/time.h>

typedef struct acados_timer_ {
    struct timeval tic;
    struct timeval toc;
} acados_timer;

#else
/* ANSI C */

/** A structure for keeping internal timer data. */
typedef struct acados_timer_ {
    struct timespec tic;
    struct timespec toc;
} acados_timer;

#endif /* __STDC_VERSION__ >= 199901L */

#endif /* (defined _WIN32 || defined _WIN64) */

/** A function for measurement of the current time. */
void acados_tic(acados_timer* t);

/** A function which returns the elapsed time. */
real_t acados_toc(acados_timer* t);

#endif

#endif  // ACADOS_UTILS_TIMING_H_

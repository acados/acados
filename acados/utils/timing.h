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

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MEASURE_TIMINGS
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

        #if __STDC_VERSION__ >= 199901L  // C99 Mode

            #include <sys/stat.h>
            #include <sys/time.h>

            typedef struct acados_timer_ {
                struct timeval tic;
                struct timeval toc;
            } acados_timer;

        #else  // ANSI C Mode

            /** A structure for keeping internal timer data. */
            typedef struct acados_timer_ {
                struct timespec tic;
                struct timespec toc;
            } acados_timer;

        #endif  // __STDC_VERSION__ >= 199901L

    #endif  // (defined _WIN32 || defined _WIN64)

#else

    // Dummy type when timings are off
    typedef real_t acados_timer;

#endif  // MEASURE_TIMINGS

/** A function for measurement of the current time. */
void acados_tic(acados_timer* t);

/** A function which returns the elapsed time. */
real_t acados_toc(acados_timer* t);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_TIMING_H_

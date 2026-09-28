/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_UTILS_TIMING_H_
#define ACADOS_UTILS_TIMING_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

#if (defined _WIN32 || defined _WIN64)

/* Use Windows QueryPerformanceCounter for timing. */
#include <Windows.h>

/** A structure for keeping internal timer data. */
typedef struct acados_timer_
{
    LARGE_INTEGER tic;
    LARGE_INTEGER toc;
    LARGE_INTEGER freq;
} acados_timer;

#elif defined(__APPLE__)

#include <mach/mach_time.h>

/** A structure for keeping internal timer data. */
typedef struct acados_timer_
{
    uint64_t tic;
    uint64_t toc;
    mach_timebase_info_data_t tinfo;
} acados_timer;

#elif defined(__MABX2__)

#include <brtenv.h>

typedef struct acados_timer_
{
    double time;
} acados_timer;

#elif defined(_DS1104)

#include <brtenv.h>

typedef struct acados_timer_
{
    double time;
} acados_timer;

#else

/* Use POSIX clock_gettime() for timing on non-Windows machines. */
#include <time.h>

#if (__STDC_VERSION__ >= 199901L) && !(defined __MINGW32__ || defined __MINGW64__)  // C99 Mode

#include <sys/stat.h>
#include <sys/time.h>

typedef struct acados_timer_
{
    struct timeval tic;
    struct timeval toc;
} acados_timer;

#else  // ANSI C Mode

/** A structure for keeping internal timer data. */
typedef struct acados_timer_
{
    struct timespec tic;
    struct timespec toc;
} acados_timer;

#endif  // __STDC_VERSION__ >= 199901L

#endif  // (defined _WIN32 || defined _WIN64)

/** A function for measurement of the current time. */
void acados_tic(acados_timer* t);

/** A function which returns the elapsed time. */
real_t acados_toc(acados_timer* t);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_TIMING_H_

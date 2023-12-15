/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include "acados/utils/timing.h"


#if (defined _WIN32 || defined _WIN64) && !(defined __MINGW32__ || defined __MINGW64__)

void acados_tic(acados_timer* t)
{
    QueryPerformanceFrequency(&t->freq);
    QueryPerformanceCounter(&t->tic);
}

real_t acados_toc(acados_timer* t)
{
    QueryPerformanceCounter(&t->toc);
    return ((t->toc.QuadPart - t->tic.QuadPart) / (real_t) t->freq.QuadPart);
}

#elif defined(__APPLE__)
void acados_tic(acados_timer* t)
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

real_t acados_toc(acados_timer* t)
{
    uint64_t duration; /* elapsed time in clock cycles*/

    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (real_t) duration / 1e9;
}

#elif defined(__MABX2__)

void acados_tic(acados_timer* t)
{
    ds1401_tic_start();
    t->time = ds1401_tic_read();
}

real_t acados_toc(acados_timer* t) { return ds1401_tic_read() - t->time; }

#else

#if (__STDC_VERSION__ >= 199901L) && !(defined __MINGW32__ || defined __MINGW64__) // C99 Mode

/* read current time */
void acados_tic(acados_timer* t) { gettimeofday(&t->tic, 0); }
/* return time passed since last call to tic on this timer */
real_t acados_toc(acados_timer* t)
{
    struct timeval temp;

    gettimeofday(&t->toc, 0);

    if ((t->toc.tv_usec - t->tic.tv_usec) < 0)
    {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
        temp.tv_usec = 1000000 + t->toc.tv_usec - t->tic.tv_usec;
    }
    else
    {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
        temp.tv_usec = t->toc.tv_usec - t->tic.tv_usec;
    }

    return (real_t) temp.tv_sec + (real_t) temp.tv_usec / 1e6;
}

#else  // ANSI C Mode

/* read current time */
void acados_tic(acados_timer* t) { clock_gettime(CLOCK_MONOTONIC, &t->tic); }
/* return time passed since last call to tic on this timer */
real_t acados_toc(acados_timer* t)
{
    struct timespec temp;

    clock_gettime(CLOCK_MONOTONIC, &t->toc);

    if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0)
    {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
        temp.tv_nsec = 1000000000 + t->toc.tv_nsec - t->tic.tv_nsec;
    }
    else
    {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
        temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
    }

    return (real_t) temp.tv_sec + (real_t) temp.tv_nsec / 1e9;
}

#endif  // __STDC_VERSION__ >= 199901L

#endif  // (defined _WIN32 || _WIN64)

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

#include "acados/utils/timing.h"

#ifdef MEASURE_TIMINGS

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

#elif defined(__DSPACE__)

void acados_tic(acados_timer* t)
{
    ds1401_tic_start();
    t->time = ds1401_tic_read();
}

real_t acados_toc(acados_timer* t) { return ds1401_tic_read() - t->time; }

#else

#if __STDC_VERSION__ >= 199901L  // C99 Mode

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

#else  // Dummy functions when timing is off

void acados_tic(acados_timer *t) {}
real_t acados_toc(acados_timer *t) { return 0; }

#endif  // MEASURE_TIMINGS

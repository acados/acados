
#include "common_header.h"
#include "timing_functions.h"

#include <stdio.h>

#if (defined WIN32 || _WIN64)

void tic(timer* t)
{
	QueryPerformanceFrequency(&t->freq);
	QueryPerformanceCounter(&t->tic);
}

real_t toc(timer* t)
{
	QueryPerformanceCounter(&t->toc);
	return ((t->toc.QuadPart - t->tic.QuadPart) / (real_t)t->freq.QuadPart);
}


#elif (defined __APPLE__)

void tic(timer* t)
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

real_t toc(timer* t)
{

    uint64_t duration; /* elapsed time in clock cycles*/

    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (real_t)duration / 1e9;
}

#else

#if __STDC_VERSION__ >= 199901L
/* C99 mode */

/* read current time */
void tic(timer* t)
{
	gettimeofday(&t->tic, 0);
}

/* return time passed since last call to tic on this timer */
real_t toc(timer* t)
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

	return (real_t)temp.tv_sec + (real_t)temp.tv_usec / 1e6;
}

#else
/* ANSI */

/* read current time */
void tic(timer* t)
{
	clock_gettime(CLOCK_MONOTONIC, &t->tic);
}


/* return time passed since last call to tic on this timer */
real_t toc(timer* t)
{
	struct timespec temp;

	clock_gettime(CLOCK_MONOTONIC, &t->toc);

	if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0)
	{
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
		temp.tv_nsec = 1000000000+t->toc.tv_nsec - t->tic.tv_nsec;
	}
	else
	{
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
		temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
	}

	return (real_t)temp.tv_sec + (real_t)temp.tv_nsec / 1e9;
}

#endif /* __STDC_VERSION__ >= 199901L */

#endif /* (defined WIN32 || _WIN64) */

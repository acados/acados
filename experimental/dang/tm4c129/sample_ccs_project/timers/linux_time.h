/*
 * linux_time.h
 *
 *  Created on: May 2, 2017
 *      Author: Dang Doan
 */

#ifndef LINUX_TIME_H_
#define LINUX_TIME_H_

#include <stdint.h>
#include "clock-arch.h"

//
// Define struct type in Linux style
//
typedef struct timeval {
    long    tv_sec;     /* seconds */
    long    tv_usec;    /* and microseconds */
} timeval;

//
// The interrupt handler for the SysTick interrupt
//
void SysTickIntHandler(void);

//
// Function that mimics gettimeofday in Linux kernel
//
int gettimeofday(struct timeval *, void *);

#endif /* LINUX_TIME_H_ */

/*
 * linux_time.c
 *
 *  Created on: May 2, 2017
 *      Author: Dang Doan

 * This function is a replacement for gettimeofday() in the Linux setting
 * and will not provide absolute time.
 * The use of this function is ONLY meaningful in comparing time, like using
 * acado_tic() and acado_toc() to measure the elapsed time.

 * Note: this function uses the SysTick feature of ARM Cortex-M4.
 * The SysTick and the SysTick interrupt handler need to be enabled.
 * In main program, SysTick period needs to be set such that 
 * there are CLOCK_CONF_SECOND SysTick cycles (interrupts) per second.
 * i.e. SysTick timer counts down from full to 0 in 1/CLOCK_CONF_SECOND s.
 * The function uses g_ui32SysClock as the system clock rate
 * (for dev. board TM4C129: 120 MHz),
 * and g_ui32TickCounter as the number of interrupt occurrences.

 * Put this code in the main program:
 *  + Global variable:
 //****************************************************************************
 //
 // System clock rate in Hz.
 //
 uint32_t g_ui32SysClock;
 //****************************************************************************
 *  + In function main():
 //****************************************************************************
    //
    // Set the clocking to run directly from the crystal at 120MHz.
    //
    g_ui32SysClock = MAP_SysCtlClockFreqSet((SYSCTL_XTAL_25MHZ |
                                             SYSCTL_OSC_MAIN |
                                             SYSCTL_USE_PLL |
                                             SYSCTL_CFG_VCO_480), 120000000);

    SysTickPeriodSet(g_ui32SysClock/CLOCK_CONF_SECOND-1);

    //
    // Enable interrupts to the processor.
    //
    IntMasterEnable();

    //
    // Enable the SysTick Interrupt.
    //
    SysTickIntEnable();

    //
    // Enable SysTick.
    //
    SysTickEnable();
 //****************************************************************************
*/

#include "linux_time.h"

extern uint32_t g_ui32SysClock;
extern uint32_t g_ui32TickCounter;
extern uint32_t SysTickValueGet();

void SysTickIntHandler(void)
{
    //
    // Increment the system tick count.
    //
    g_ui32TickCounter++;
    //
    // Display SysTick counter every second (for debug)
    //
    //if (g_ui32TickCounter % CLOCK_CONF_SECOND == 0)
    //{
    //    UARTprintf("%d\n",g_ui32TickCounter);
    //}
}

//
// Function that mimics gettimeofday in Linux kernel
//
int gettimeofday(struct timeval *tv, void *tz) {
    uint32_t ui32TempSysTickValue;
    uint32_t ui32CyclesPerMicrosec;
    uint32_t ui32usecPerRound = 1000000/CLOCK_CONF_SECOND;
    ui32TempSysTickValue = SysTickValueGet();
    // TO DO: Better stop the counter before reading SysTick value and SysTick counter
    //        to avoid counting at the transition 0 -> full, that could cause 1 SysTickcounter inaccuracy
    tv->tv_sec = g_ui32TickCounter / CLOCK_CONF_SECOND;
    ui32CyclesPerMicrosec = g_ui32SysClock / 1000000;
    // ui32TempSysTickValue / ui32CyclesPerMicrosec = number of microseconds to count down to finish that SysTick round
    // One SysTick round = 1/CLOCK_CONF_SECOND sec = 1000000/CLOCK_CONF_SECOND usec
    // Number of usec/round - #usec_to_finish_the_round = #usec_elapsed_in_that_round
    tv->tv_usec = (g_ui32TickCounter % CLOCK_CONF_SECOND)*ui32usecPerRound + (ui32usecPerRound - ui32TempSysTickValue / ui32CyclesPerMicrosec);
    return 0;
}

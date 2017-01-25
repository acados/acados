#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wimplicit-function-declaration"

/*----------------------------------------------------------------------------
                                    E.S.O.
 -----------------------------------------------------------------------------
   File name    :   gnuplot_i.h
   Author       :   N. Devillard
   Created on   :   Fri Sept 26 1997
   Software     :   ANSI C under Solaris Unix
                    Part of ECLIPSE library for Adonis
   Description  :   C interface to gnuplot

    gnuplot is a freely available, command-driven graphical display tool for
    Unix. It compiles and works quite well on a number of Unix flavours as
    well as other operating systems. The following module enables sending
    display requests to gnuplot through simple C calls.

 ---------------------------------------------------------------------------*/

/*

 $Id: gnuplot_i.h,v 1.2 2000/04/18 12:31:17 ndevilla Exp $
 $Author: ndevilla $
 $Date: 2000/04/18 12:31:17 $
 $Revision: 1.2 $

 */

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>

#define GP_MAX_TMP_FILES    64
#define GP_TMP_NAME_SIZE    512
#define GP_CMD_SIZE     	1024

/*---------------------------------------------------------------------------
                                New Types
 ---------------------------------------------------------------------------*/
/*
 * This structure holds all necessary information to talk to a gnuplot
 * session.
 */
typedef struct _GNUPLOT_CTRL_ {
    /* command file handling */
    FILE    * gnucmd ;

    /* Plotting options */
    int       nplots ;      /* Number of active plots at the moment */
    char      pstyle[32] ;  /* Current plotting style */

    /* temporary files opened */
    char      to_delete[GP_MAX_TMP_FILES][GP_TMP_NAME_SIZE] ;
    int       ntmp ;


} gnuplot_ctrl ;


#ifndef _ECLIPSE_TYPES_H_
/*
 * dpoint is convenient to store signals which have definition both on x and
 * y axis.
 */
typedef struct _DPOINT_ {
    double  x ;
    double  y ;
} dpoint ;
#endif


/*---------------------------------------------------------------------------
                        Function ANSI C prototypes
 ---------------------------------------------------------------------------*/
int check_X_display(int activate);
char * gnuplot_get_program_path(char * pname);
gnuplot_ctrl * gnuplot_init(void);
void gnuplot_close(gnuplot_ctrl * handle);
void gnuplot_cmd ( gnuplot_ctrl *  handle, char *  cmd, ...);
void gnuplot_setstyle(gnuplot_ctrl * h, char * plot_style);
void gnuplot_set_xlabel ( gnuplot_ctrl * h, char * label );
void gnuplot_set_ylabel(gnuplot_ctrl * h, char * label);
void gnuplot_resetplot(gnuplot_ctrl * h);
void gnuplot_plot1d_var1(gnuplot_ctrl*handle,double *d, int n_point, char *title ) ;
void gnuplot_plot1d_var2(gnuplot_ctrl *handle, dpoint *d, int n_points, char *title );
void gnuplot_plot1d_var2v (gnuplot_ctrl *handle, double *x, double *y,
    int n_points, char *title );
void gnuplot_plot_slope(gnuplot_ctrl *handle, double a, double b, char *title ) ;
void gnuplot_plot_equation(gnuplot_ctrl *h, char *equation, char *title ) ;

#endif

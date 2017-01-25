#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wint-conversion"

/*----------------------------------------------------------------------------
                                    E.S.O.
 -----------------------------------------------------------------------------
   File name    :   gnuplot_i.c
   Author       :   N. Devillard
   Created on   :   Fri Sept 26 1997
   Language     :   ANSI C
   Description  :   C interface to gnuplot

    gnuplot is a freely available, command-driven graphical display tool for
    Unix. It compiles and works quite well on a number of Unix flavours as
    well as other operating systems. The following module enables sending
    display requests to gnuplot through simple C calls.

 ---------------------------------------------------------------------------*/

/*

 $Id: gnuplot_i.c,v 1.3 2000/05/11 11:00:57 ndevilla Exp $
 $Author: ndevilla $
 $Date: 2000/05/11 11:00:57 $
 $Revision: 1.3 $

 */

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/

#include "gnuplot_i.h"

/*---------------------------------------------------------------------------
                                Defines
 ---------------------------------------------------------------------------*/

/* Maximal size of a gnuplot command */
#define GP_CMD_SIZE     1024
/* Maximal size of a plot title */
#define GP_TITLE_SIZE   80
/* Maximal size for an equation */
#define GP_EQ_SIZE      512


/*---------------------------------------------------------------------------
                            Function codes
 ---------------------------------------------------------------------------*/




int check_X_display ( int activate )
/*-------------------------------------------------------------------------*/
/*
  @name		check_X_display
  @memo		Checks out if the DISPLAY environment variable is set.
  @param	activate int flag
  @return	int 1 if the variable is set, 0 otherwise.
  @doc

  This function checks out the DISPLAY environment variable to see if
  it exists. It does not check if the display is actually correctly
  configured. If you do not want to activate this check (e.g. on
  systems that do not support this kind of display mechanism), pass a
  0 integer as the activate flag. Any other value will activate it.
*/
/*--------------------------------------------------------------------------*/
{
  char *display;

  if ( !activate )
  {
    return 1;
  }

  display = getenv ( "DISPLAY" );

  if ( display == NULL )
  {
    fprintf ( stderr, "cannot find DISPLAY variable: is it set?\n");
    return 1 ;
  }
  else
  {
    return 0 ;
  }
}

/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_get_program_path
  @memo		Find out where a command lives in your PATH.
  @param	pname Name of the program to look for.
  @return	pointer to statically allocated character string.
  @doc

  This is the C equivalent to the 'which' command in Unix. It parses
  out your PATH environment variable to find out where a command
  lives. The returned character string is statically allocated within
  this function, i.e. there is no need to free it. Beware that the
  contents of this string will change from one call to the next,
  though (as all static variables in a function).

  The input character string must be the name of a command without
  prefixing path of any kind, i.e. only the command name. The returned
  string is the path in which a command matching the same name was
  found.

  Examples (assuming there is a prog named 'hello' in the cwd):

  \begin{itemize}
  \item gnuplot_get_program_path("hello") returns "."
  \item gnuplot_get_program_path("ls") returns "/bin"
  \item gnuplot_get_program_path("csh") returns "/usr/bin"
  \item gnuplot_get_program_path("/bin/ls") returns NULL
  \end{itemize}

 */
/*-------------------------------------------------------------------------*/

#define MAXNAMESZ       4096
char * gnuplot_get_program_path(char * pname)
{
    int         i, j, lg;
    char    *   path;
    static char buf[MAXNAMESZ];

    /* Trivial case: try in CWD */
    sprintf(buf, "./%s", pname) ;
    if (access(buf, X_OK)==0) {
        sprintf(buf, ".");
        return buf ;
    }
    /* Try out in all paths given in the PATH variable */
    buf[0] = 0;
    path = getenv("PATH") ;
    if (path!=NULL) {
        for (i=0; path[i]; ) {
            for (j=i ; (path[j]) && (path[j]!=':') ; j++);
            lg = j - i;
            strncpy(buf, path + i, lg);
            if (lg == 0) buf[lg++] = '.';
            buf[lg++] = '/';
            strcpy(buf + lg, pname);
            if (access(buf, X_OK) == 0) {
                /* Found it! */
                break ;
            }
            buf[0] = 0;
            i = j;
            if (path[i] == ':') i++ ;
        }
    }
    /* If the buffer is still empty, the command was not found */
    if (buf[0] == 0) return NULL ;
    /* Otherwise truncate the command name to yield path only */
    lg = strlen(buf) - 1 ;
    while (buf[lg]!='/') {
        buf[lg]=0 ;
        lg -- ;
    }
    buf[lg] = 0;
    return buf ;
}
#undef MAXNAMESZ



/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_init
  @memo		Opens up a gnuplot session, ready to receive commands.
  @return	Newly allocated gnuplot control structure.
  @doc

  This opens up a new gnuplot session, ready for input. The struct
  controlling a gnuplot session should remain opaque and only be
  accessed through the provided functions.
 */
/*--------------------------------------------------------------------------*/

gnuplot_ctrl * gnuplot_init ( void )
{
    gnuplot_ctrl *  handle ;

    if (check_X_display(1)) return NULL ;

	if (gnuplot_get_program_path("gnuplot")==NULL) {
	  fprintf(stderr, "cannot find gnuplot in your PATH");
	  return NULL ;
	}

    /*
     * Structure initialization:
     */
    handle = ( gnuplot_ctrl * ) malloc ( sizeof ( gnuplot_ctrl ) ) ;
    handle->nplots = 0 ;
    gnuplot_setstyle(handle, "points") ;
    handle->ntmp = 0 ;

    handle->gnucmd = popen("gnuplot", "w") ;
    if (handle->gnucmd == NULL) {
        fprintf(stderr, "error starting gnuplot\n") ;
        free(handle) ;
        return NULL ;
    }
    return handle;
}


/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_close
  @memo		Closes a gnuplot session previously opened by gnuplot_init()
  @param	handle Gnuplot session control handle.
  @return	void
  @doc

  Kills the child PID and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.

 */
/*--------------------------------------------------------------------------*/

void gnuplot_close(gnuplot_ctrl * handle)
{
    int     i ;
    if (check_X_display(1)) return ;
    if (handle->ntmp) {
        for (i=0 ; i<handle->ntmp ; i++) {
            remove(handle->to_delete[i]) ;
        }
    }
    if (pclose(handle->gnucmd) == -1) {
        fprintf(stderr, "problem closing communication to gnuplot\n") ;
        return ;
    }
    free(handle) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_cmd
  @memo		Sends a command to an active gnuplot session.
  @param	handle Gnuplot session control handle
  @param	cmd    Command to send, same as a printf statement.
  @return	void
  @doc

  This sends a string to an active gnuplot session, to be executed.
  There is strictly no way to know if the command has been
  successfully executed or not.
  The command syntax is the same as printf.

  Examples:

  \begin{itemize}
  \item gnuplot_cmd(g, "plot %d*x", 23.0);
  \item gnuplot_cmd(g, "plot %g * cos(%g * x)", 32.0, -3.0);
  \end{itemize}

 */
/*--------------------------------------------------------------------------*/

void gnuplot_cmd ( gnuplot_ctrl *  handle, char *  cmd, ...)
{
    va_list ap ;
    char    local_cmd[GP_CMD_SIZE];

    va_start(ap, cmd);
    vsprintf(local_cmd, cmd, ap);
    va_end(ap);

    strcat(local_cmd, "\n");

    fputs(local_cmd, handle->gnucmd) ;
    fflush(handle->gnucmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_setstyle
  @memo		Change the plotting style of a gnuplot session.
  @param	h Gnuplot session control handle
  @param	plot_style Plotting-style to use (character string)
  @return	void
  @doc

  The provided plotting style is a character string. It must be one of
  the following:

  \begin{itemize}
  \item {\it lines}
  \item {\it points}
  \item {\it linespoints}
  \item {\it impulses}
  \item {\it dots}
  \item {\it steps}
  \item {\it errorbars}
  \item {\it boxes}
  \item {\it boxeserrorbars}
  \end{itemize}
 */
/*--------------------------------------------------------------------------*/

void gnuplot_setstyle(gnuplot_ctrl * h, char * plot_style)
{
    if (strcmp(plot_style, "lines") &&
        strcmp(plot_style, "points") &&
        strcmp(plot_style, "linespoints") &&
        strcmp(plot_style, "impulses") &&
        strcmp(plot_style, "dots") &&
        strcmp(plot_style, "steps") &&
        strcmp(plot_style, "errorbars") &&
        strcmp(plot_style, "boxes") &&
        strcmp(plot_style, "boxerrorbars")) {
        fprintf(stderr, "warning: requested style unknown: therefore, using points\n") ;
        (void)strcpy(h->pstyle, "points") ;
    } else {
        (void)strcpy(h->pstyle, plot_style) ;
    }
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_set_xlabel
  @memo		Sets the x label of a gnuplot session.
  @param	h Gnuplot session control handle.
  @param	label Character string to use for X label.
  @return	void
  @doc

  Sets the x label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_xlabel ( gnuplot_ctrl * h, char * label )
{
    char    cmd[GP_CMD_SIZE] ;

    (void)sprintf(cmd, "set xlabel \"%s\"", label) ;
    gnuplot_cmd(h, cmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_set_ylabel
  @memo		Sets the y label of a gnuplot session.
  @param	h Gnuplot session control handle.
  @param	label Character string to use for Y label.
  @return	void
  @doc

  Sets the y label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_ylabel(gnuplot_ctrl * h, char * label)
{
    char    cmd[GP_CMD_SIZE] ;

    (void)sprintf(cmd, "set ylabel \"%s\"", label) ;
    gnuplot_cmd(h, cmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_resetplot
  @memo		Resets a gnuplot session (next plot will erase previous ones).
  @param	h Gnuplot session control handle.
  @return	void
  @doc

  Resets a gnuplot session, i.e. the next plot will erase all previous
  ones.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_resetplot(gnuplot_ctrl * h)
{
    int     i ;
    if (h->ntmp) {
        for (i=0 ; i<h->ntmp ; i++) {
            remove(h->to_delete[i]) ;
        }
    }
    h->ntmp = 0 ;
    h->nplots = 0 ;
    return ;
}



/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_plot1d_var1
  @memo		Plots a 2d graph from a list of doubles.
  @param	handle	Gnuplot session control handle.
  @param	d		Pointer to a list of doubles.
  @param	n_point	Number of doubles in the list.
  @param	title	Title of the plot.
  @return	void
  @doc

  Plots out a 2d graph from a list of doubles. The x-coordinate is the
  index of the double in the list, the y coordinate is the double in
  the list.

  Example:

  \begin{verbatim}
    gnuplot_ctrl    *h ;
    double          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i] = (double)(i*i) ;
    }
    gnuplot_plot1d_var1(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  \end{verbatim}

 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot1d_var1(
    gnuplot_ctrl    *   handle,
    double          *   d,
    int                 n_point,
    char            *   title
)
{
    int         i ;
    FILE    *   tmp ;
    char    *   name ;
    char        cmd[GP_CMD_SIZE] ;
    char        line[GP_CMD_SIZE] ;

    /* can we open one more temporary file? */
    if (handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return ;
    }

    /* Open temporary file for output   */
    if ((name = tmpnam(NULL)) == (char*)NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }
    if ((tmp = fopen(name, "w")) == NULL) {
        fprintf(stderr, "cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Store file name in array for future deletion */
    (void)strcpy(handle->to_delete[handle->ntmp], name) ;
    handle->ntmp ++ ;

    /* Write data to this file  */
    for (i=0 ; i<n_point ; i++) {
        (void)fprintf(tmp, "%g\n", d[i]) ;
    }
    (void)fflush(tmp) ;
    (void)fclose(tmp) ;

    /* Command to be sent to gnuplot    */
    if (handle->nplots > 0) {
        (void)strcpy(cmd, "replot") ;
    } else {
        (void)strcpy(cmd, "plot") ;
    }

    if (title == NULL) {
        (void)sprintf(line, "%s \"%s\" with %s", cmd, name, handle->pstyle) ;
    } else {
        (void)sprintf(line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, handle->pstyle) ;
    }

    /* send command to gnuplot  */
    gnuplot_cmd(handle, line) ;
    handle->nplots++ ;
    return ;
}



/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_plot1d_var2
  @memo		Plot a 2d graph from a list of dpoint.
  @param	handle		Gnuplot session control handle.
  @param	d			Pointer to a list of doubles.
  @param	n_points	Number of doubles in the list.
  @param	title		Title of the plot.
  @return	void
  @doc

  Plots out a 2d graph from a list of dpoints. A dpoint is a struct
  containing two fields x and y (doubles) which are plotted as they
  are on the gnuplot session.

  \begin{verbatim}
    gnuplot_ctrl    *h ;
    dpoint          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i].x = (double)(i)/10.0 ;
        d[i].y = d[i].x * d[i].x ;
    }
    gnuplot_plot1d_var2(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  \end{verbatim}
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot1d_var2(
    gnuplot_ctrl    *   handle,
    dpoint          *   d,
    int                 n_points,
    char            *   title
)
{
    int         i ;
    FILE    *   tmp ;
    char    *   name ;
    char        cmd[GP_CMD_SIZE] ;
    char        line[GP_CMD_SIZE] ;

    /* can we open one more temporary file? */
    if (handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return ;
    }

    /* Open temporary file for output   */
    if ((name = tmpnam(NULL)) == (char*)NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }
    if ((tmp = fopen(name, "w")) == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Store file name in array for future deletion */
    (void)strcpy(handle->to_delete[handle->ntmp], name) ;
    handle->ntmp ++ ;

    /* Write data to this file  */
    for (i=0 ; i<n_points ; i++) {
        (void)fprintf(tmp, "%g %g\n", d[i].x, d[i].y) ;
    }
    (void)fflush(tmp) ;
    (void)fclose(tmp) ;

    /* Command to be sent to gnuplot    */
    if (handle->nplots > 0) {
        (void)strcpy(cmd, "replot") ;
    } else {
        (void)strcpy(cmd, "plot") ;
    }

    if (title == NULL) {
        (void)sprintf(line, "%s \"%s\" with %s", cmd, name, handle->pstyle) ;
    } else {
        (void)sprintf(line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, handle->pstyle) ;
    }

    /* send command to gnuplot  */
    gnuplot_cmd(handle, line) ;
    handle->nplots++ ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_plot1d_var2v
  @memo		Plot a 2d graph from a list of pairs of points.
  @param	handle		Gnuplot session control handle.
  @param	x			The vector of x values;
  @param    y           The vector of y values;
  @param	n_points	Number of doubles in the list.
  @param	title		Title of the plot.
  @return	void
  @doc

  Plots out a 2d graph from a list of pairs of points.

  \begin{verbatim}
    gnuplot_ctrl    *h ;
    double          x[50];
    double          y[50];
    int             j;

    h = gnuplot_init() ;
    for (j=0 ; j<50 ; j++) {
        x[j] = (double)(i)/10.0 ;
        y[j] = x[j] * x[j];
    }
    gnuplot_plot1d_var2v(h, x, y, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  \end{verbatim}
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot1d_var2v (
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    int                 n_points,
    char            *   title
)
{
  int j;
  FILE *tmp;
  char *name;
  char cmd[GP_CMD_SIZE];
  char line[GP_CMD_SIZE];
/*
  Can we open one more temporary file?
*/
  if ( handle->ntmp == GP_MAX_TMP_FILES - 1 )
  {
    fprintf ( stderr,
      "Maximum # of temporary files reached (%d): cannot open more,",
      GP_MAX_TMP_FILES );
    return;
  }
/*
  Open temporary file for output
*/
  if ( ( name = tmpnam(NULL) ) == (char*) NULL )
  {
    fprintf ( stderr, "cannot create temporary file: exiting plot" );
    return ;
  }

  if ( ( tmp = fopen ( name, "w" ) ) == NULL )
  {
    fprintf ( stderr, "cannot create temporary file: exiting plot" );
    return ;
  }
/*
  Store file name in array for future deletion
*/
  strcpy ( handle->to_delete[handle->ntmp], name );
  handle->ntmp++;
/*
  Write data to this file
*/
  for ( j = 0; j < n_points; j++ )
  {
    fprintf ( tmp, "%g %g\n", x[j], y[j] ) ;
  }
  fflush ( tmp );
  fclose ( tmp );
/*
  Command to be sent to gnuplot.
*/
  if ( handle->nplots > 0 )
  {
    strcpy ( cmd, "replot" );
  }
  else
  {
    strcpy ( cmd, "plot" );
  }

  if ( title == NULL )
  {
    sprintf ( line, "%s \"%s\" with %s", cmd, name, handle->pstyle );
  }
  else
  {
    sprintf ( line, "%s \"%s\" title \"%s\" with %s", cmd, name,
                      title, handle->pstyle) ;
  }
/*
  Send command to gnuplot
*/
  gnuplot_cmd ( handle, line );
  handle->nplots++;
  return;
}
/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_plot_slope
  @memo		Plot a slope on a gnuplot session.
  @param	handle		Gnuplot session control handle.
  @param	a			Slope.
  @param	b			Intercept.
  @param	title		Title of the plot.
  @return	void
  @doc

  Plot a slope on a gnuplot session. The provided slope has an
  equation of the form:

  \begin{verbatim}
  y = ax+b
  \end{verbatim}

  Example:

  \begin{verbatim}
    gnuplot_ctrl    *   h ;
    double              a, b ;

    h = gnuplot_init() ;
    gnuplot_plot_slope(h, 1.0, 0.0, "unity slope") ;
    sleep(2) ;
    gnuplot_close(h) ;
  \end{verbatim}

 */
/*--------------------------------------------------------------------------*/


void gnuplot_plot_slope(
    gnuplot_ctrl    *   handle,
    double              a,
    double              b,
    char            *   title
)
{
    char    stitle[GP_TITLE_SIZE] ;
    char    cmd[GP_CMD_SIZE] ;

    if (title == NULL) {
        (void)strcpy(stitle, "no title") ;
    } else {
        (void)strcpy(stitle, title) ;
    }

    if (handle->nplots > 0) {
        (void)sprintf(cmd, "replot %g * x + %g title \"%s\" with %s",
                      a, b, title, handle->pstyle) ;
    } else {
        (void)sprintf(cmd, "plot %g * x + %g title \"%s\" with %s",
                      a, b, title, handle->pstyle) ;
    }
    gnuplot_cmd(handle, cmd) ;
    handle->nplots++ ;
    return ;
}



void gnuplot_plot_equation ( gnuplot_ctrl *h, char *equation, char *title )
/*-------------------------------------------------------------------------*/
/**
  @name		gnuplot_plot_equation
  @memo		Plot a curve of given equation y=f(x).
  @param	h			Gnuplot session control handle.
  @param	equation	Equation to plot.
  @param	title		Title of the plot.
  @return	void
  @doc

  Plots out a curve of given equation. The general form of the
  equation is y=f(x), you only provide the f(x) side of the equation.

  Example:

  \begin{verbatim}
        gnuplot_ctrl    *h ;
        char            eq[80] ;

        h = gnuplot_init() ;
        strcpy(eq, "sin(x) * cos(2*x)") ;
        gnuplot_plot_equation(h, eq, "sine wave", normal) ;
        gnuplot_close(h) ;
  \end{verbatim}

 */
/*--------------------------------------------------------------------------*/

{
  char  cmd[GP_CMD_SIZE];
  char  plot_str[GP_EQ_SIZE];
  char  title_str[GP_TITLE_SIZE];

  if ( title == NULL )
  {
    ( void ) strcpy ( title_str, "no title" );
  }
  else
  {
    ( void ) strcpy ( title_str, title);
  }

  if ( h->nplots > 0 )
  {
    ( void ) strcpy ( plot_str, "replot" );
  }
  else
  {
    ( void ) strcpy ( plot_str, "plot" );
  }

  ( void ) sprintf ( cmd, "%s %s title \"%s\" with %s",
    plot_str, equation, title_str, h->pstyle );

  gnuplot_cmd ( h, cmd );

  h->nplots++;

  return;
}

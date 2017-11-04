/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file src/Utils.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches (thanks to Eckhard Arnold)
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Implementation of some utilities for working with the different QProblem classes.
 */


#if defined(__WIN32__) || defined(WIN32)
  #include <windows.h>
#elif defined(LINUX) || defined(__LINUX__)
  #include <sys/stat.h>
  #include <sys/time.h>
#endif

#ifdef __MATLAB__
  #include "mex.h"
#endif


#include <qpOASES_e/Utils.h>


#ifdef __NO_SNPRINTF__
#if (!defined(_MSC_VER)) || defined(__DSPACE__) || defined(__XPCTARGET__)
/* If snprintf is not available, provide an empty implementation. */
int snprintf( char* s, size_t n, const char* format, ... )
{
	if ( n > 0 )
		s[0] = '\0';

	return 0;
}
#endif
#endif /* __NO_SNPRINTF__ */


BEGIN_NAMESPACE_QPOASES


/*
 *	p r i n t
 */
returnValue qpOASES_printV(	const real_t* const v,
							int n
							)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int i;
	char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	/* Print a vector. */
	for( i=0; i<n; ++i )
	{
		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH," %.16e\t", v[i] );
		qpOASES_myPrintf( myPrintfString );
	}
	qpOASES_myPrintf( "\n" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue qpOASES_printPV(	const real_t* const v, int n,
								const int* const V_idx
								)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int i;
	char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	/* Print a permuted vector. */
	qpOASES_myPrintf( "\t" );
	for( i=0; i<n; ++i )
	{
		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH," %.16e\t", v[ V_idx[i] ] );
		qpOASES_myPrintf( myPrintfString );
	}
	qpOASES_myPrintf( "\n" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue qpOASES_printNV(	const real_t* const v, int n,
								const char* name
								)
{
	#ifndef __SUPPRESSANYOUTPUT__

	char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	/* Print vector name ... */
	if ( name != 0 )
	{
		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"%s = ", name );
		qpOASES_myPrintf( myPrintfString );
	}

	#endif /* __SUPPRESSANYOUTPUT__ */

	/* ... and the vector itself. */
	return qpOASES_printV( v, n );
}

/*
 *	p r i n t
 */
returnValue qpOASES_printM(	const real_t* const M,
							int nrow,
							int ncol
							)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int i;

	/* Print a matrix as a collection of row vectors. */
	for( i=0; i<nrow; ++i )
		qpOASES_printV( &(M[i*ncol]), ncol );
	qpOASES_myPrintf( "\n" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue qpOASES_printPM(	const real_t* const M, int nrow, int ncol,
								const int* const ROW_idx, const int* const COL_idx
								)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int i;

	/* Print a permuted matrix as a collection of permuted row vectors. */
	for( i=0; i<nrow; ++i )
		qpOASES_printPV( &( M[ ROW_idx[i]*ncol ] ), ncol, COL_idx );
	qpOASES_myPrintf( "\n" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue qpOASES_printNM(	const real_t* const M, int nrow, int ncol,
								const char* name
								)
{
	#ifndef __SUPPRESSANYOUTPUT__

	char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	/* Print matrix name ... */
	if ( name != 0 )
	{
		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"%s = ", name );
		qpOASES_myPrintf( myPrintfString );
	}

	#endif /* __SUPPRESSANYOUTPUT__ */

	/* ... and the matrix itself. */
	return qpOASES_printM( M, nrow, ncol );
}


/*
 *	p r i n t
 */
returnValue qpOASES_printI(	const int* const _index,
							int n
							)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int i;
	char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	/* Print a indexlist. */
	qpOASES_myPrintf( "\t" );
	for( i=0; i<n; ++i )
	{
		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH," %d\t", _index[i] );
		qpOASES_myPrintf( myPrintfString );
	}
	qpOASES_myPrintf( "\n" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue qpOASES_printNI(	const int* const _index, int n,
								const char* name
								)
{
	#ifndef __SUPPRESSANYOUTPUT__

	char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	/* Print indexlist name ... */
	if ( name != 0 )
	{
		snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"%s = ", name );
		qpOASES_myPrintf( myPrintfString );
	}

	#endif /* __SUPPRESSANYOUTPUT__ */

	/* ... and the indexlist itself. */
	return qpOASES_printI( _index, n );
}


/*
 *	m y P r i n t f
 */
returnValue qpOASES_myPrintf( const char* s )
{
	#ifndef __SUPPRESSANYOUTPUT__

	  #ifdef __MATLAB__
	  if ( s == 0 )
	      return RET_INVALID_ARGUMENTS;
	  mexPrintf( s );
	  #else
	  FILE* outputfile;
	  if ( s == 0 )
	      return RET_INVALID_ARGUMENTS;
	  outputfile = MessageHandling_getOutputFile( qpOASES_getGlobalMessageHandler() );
	  if ( outputfile == 0 )
		  return THROWERROR( RET_NO_GLOBAL_MESSAGE_OUTPUTFILE );
	  fprintf( outputfile, "%s", s );
	  #endif /* __MATLAB__ */

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t C o p y r i g h t N o t i c e
 */
returnValue qpOASES_printCopyrightNotice( )
{
	#ifndef __SUPPRESSANYOUTPUT__
		#ifndef __XPCTARGET__
		#ifndef __DSPACE__
		#ifndef __NO_COPYRIGHT__
		qpOASES_myPrintf( "\nqpOASES_embedded -- An Embedded Implementation of the Online Active Set Strategy.\nCopyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,\nChristian Kirches et al. All rights reserved.\n\nqpOASES is distributed under the terms of the \nGNU Lesser General Public License 2.1 in the hope that it will be \nuseful, but WITHOUT ANY WARRANTY; without even the implied warranty \nof MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. \nSee the GNU Lesser General Public License for more details.\n\n" );
		#endif /* __NO_COPYRIGHT__ */
		#endif /* __DSPACE__ */
		#endif /* __XPCTARGET__ */
	#endif /* __SUPPRESSANYOUTPUT__ */
	return SUCCESSFUL_RETURN;
}


/*
 *	r e a d F r o m F i l e
 */
returnValue qpOASES_readFromFileM(	real_t* data, int nrow, int ncol,
									const char* datafilename
									)
{
	// #ifndef __SUPPRESSANYOUTPUT__

	int i, j;
	real_t float_data;
	FILE* datafile;
	char errstr[QPOASES_MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
		return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_OPEN_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
	}

	/* 2) Read data from file. */
	for( i=0; i<nrow; ++i )
	{
		for( j=0; j<ncol; ++j )
		{
			#ifdef __USE_SINGLE_PRECISION__
			if ( fscanf( datafile, "%f ", &float_data ) == 0 )
			#else
			if ( fscanf( datafile, "%lf ", &float_data ) == 0 )
			#endif /* __USE_SINGLE_PRECISION__ */
			{
				fclose( datafile );
				snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
				return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_READ_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
			}
			data[i*ncol + j] = ( (real_t) float_data );
		}
	}

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;

	// #else

	// return RET_NOT_YET_IMPLEMENTED;

	// #endif /* __SUPPRESSANYOUTPUT__ */
}


/*
 *	r e a d F r o m F i l e
 */
returnValue qpOASES_readFromFileV(	real_t* data, int n,
									const char* datafilename
									)
{
	return qpOASES_readFromFileM( data, n, 1, datafilename );
}



/*
 *	r e a d F r o m F i l e
 */
returnValue qpOASES_readFromFileI(	int* data, int n,
									const char* datafilename
									)
{
	// #ifndef __SUPPRESSANYOUTPUT__

	int i;
	FILE* datafile;
	char errstr[QPOASES_MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
		return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_OPEN_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
	}

	/* 2) Read data from file. */
	for( i=0; i<n; ++i )
	{
		if ( fscanf( datafile, "%d\n", &(data[i]) ) == 0 )
		{
			fclose( datafile );
			snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
			return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_READ_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;

    // #else

	// return RET_NOT_YET_IMPLEMENTED;

	// #endif /* __SUPPRESSANYOUTPUT__ */
}


/*
 *	w r i t e I n t o F i l e
 */
returnValue qpOASES_writeIntoFileM(	const real_t* const data, int nrow, int ncol,
									const char* datafilename, BooleanType append
									)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int i, j;
	FILE* datafile;
	char errstr[QPOASES_MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( append == BT_TRUE )
	{
		/* append data */
		if ( ( datafile = fopen( datafilename, "a" ) ) == 0 )
		{
			snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
			return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_OPEN_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}
	else
	{
		/* do not append data */
		if ( ( datafile = fopen( datafilename, "w" ) ) == 0 )
		{
			snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
			return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_OPEN_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}

	/* 2) Write data into file. */
	for( i=0; i<nrow; ++i )
	{
		for( j=0; j<ncol; ++j )
		 	fprintf( datafile, "%.16e ", data[i*ncol+j] );

		fprintf( datafile, "\n" );
	}

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;

	#else /* __SUPPRESSANYOUTPUT__ */

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	w r i t e I n t o F i l e
 */
returnValue qpOASES_writeIntoFileV(	const real_t* const data, int n,
									const char* datafilename, BooleanType append
									)
{
	return qpOASES_writeIntoFileM( data,1,n,datafilename,append );
}


/*
 *	w r i t e I n t o F i l e
 */
returnValue qpOASES_writeIntoFileI(	const int* const integer, int n,
									const char* datafilename, BooleanType append
									)
{
	#ifndef __SUPPRESSANYOUTPUT__

	int i;

	FILE* datafile;
	char errstr[QPOASES_MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( append == BT_TRUE )
	{
		/* append data */
		if ( ( datafile = fopen( datafilename, "a" ) ) == 0 )
		{
			snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
			return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_OPEN_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}
	else
	{
		/* do not append data */
		if ( ( datafile = fopen( datafilename, "w" ) ) == 0 )
		{
			snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
			return MessageHandling_throwError( qpOASES_getGlobalMessageHandler(),RET_UNABLE_TO_OPEN_FILE,errstr,__FUNC__,__FILE__,__LINE__,VS_VISIBLE );
		}
	}

	/* 2) Write data into file. */
	for( i=0; i<n; ++i )
		fprintf( datafile, "%d\n", integer[i] );

	/* 3) Close file. */
	fclose( datafile );

	return SUCCESSFUL_RETURN;

    #else /* __SUPPRESSANYOUTPUT__ */

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	w r i t e I n t o M a t F i l e
 */
returnValue qpOASES_writeIntoMatFile(	FILE* const matFile,
										const double* const data, int nRows, int nCols, const char* name
										)
{
	/*  Note, this code snippet has been inspired from the document
	 *  "Matlab(R) MAT-file Format, R2013b" by MathWorks */

	#ifndef __SUPPRESSANYOUTPUT__

	int ii, jj;
	MatMatrixHeader var;

	if ( ( matFile == 0 ) || ( data == 0 ) || ( nRows < 0 ) || ( nCols < 0 ) || ( name == 0 ) )
		return RET_INVALID_ARGUMENTS;

	/* setup variable header */
	var.numericFormat = 0000;  /* IEEE Little Endian - reserved - double precision (64 bits) - numeric full matrix */
	var.nRows         = nRows; /* number of rows */
	var.nCols         = nCols; /* number of columns */
	var.imaginaryPart = 0;     /* no imaginary part */
	var.nCharName     = (long)(strlen(name))+1; /* matrix name length  */

	/* write variable header to mat file */
	if ( fwrite( &var, sizeof(MatMatrixHeader),1,  matFile ) < 1 )
		return RET_UNABLE_TO_WRITE_FILE;

	if ( fwrite( name, sizeof(char),(unsigned long)(var.nCharName), matFile ) < 1 )
		return RET_UNABLE_TO_WRITE_FILE;

	for ( ii=0; ii<nCols; ++ii )
		for ( jj=0; jj<nRows; ++jj )
			if ( fwrite( &(data[jj*nCols+ii]), sizeof(double),1, matFile ) < 1 )
				return RET_UNABLE_TO_WRITE_FILE;

	return SUCCESSFUL_RETURN;

    #else /* __SUPPRESSANYOUTPUT__ */

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	w r i t e I n t o M a t F i l e
 */
returnValue qpOASES_writeIntoMatFileI(	FILE* const matFile,
										const int* const data, int nRows, int nCols, const char* name
										)
{
	/*  Note, this code snippet has been inspired from the document
	 *  "Matlab(R) MAT-file Format, R2013b" by MathWorks */

	#ifndef __SUPPRESSANYOUTPUT__

	int ii, jj;
	MatMatrixHeader var;

	if ( ( matFile == 0 ) || ( data == 0 ) || ( nRows < 0 ) || ( nCols < 0 ) || ( name == 0 ) )
		return RET_INVALID_ARGUMENTS;

	/* setup variable header */
	var.numericFormat = 0000;  /* IEEE Little Endian - reserved - double precision (64 bits) - numeric full matrix */
	var.nRows         = nRows; /* number of rows */
	var.nCols         = nCols; /* number of columns */
	var.imaginaryPart = 0;     /* no imaginary part */
	var.nCharName     = (long)(strlen(name))+1; /* matrix name length  */

	/* write variable header to mat file */
	if ( fwrite( &var, sizeof(MatMatrixHeader),1,  matFile ) < 1 )
		return RET_UNABLE_TO_WRITE_FILE;

	if ( fwrite( name, sizeof(char),(unsigned long)(var.nCharName), matFile ) < 1 )
		return RET_UNABLE_TO_WRITE_FILE;

	for ( ii=0; ii<nCols; ++ii ) {
		for ( jj=0; jj<nRows; ++jj ) {
			double real_data = (double) data[jj*nCols+ii];
			if ( fwrite( &real_data, sizeof(double),1, matFile ) < 1 )
							return RET_UNABLE_TO_WRITE_FILE;
		}
	}

	return SUCCESSFUL_RETURN;

    #else /* __SUPPRESSANYOUTPUT__ */

	return RET_NOT_YET_IMPLEMENTED;

	#endif
}


/*
 *	g e t C P U t i m e
 */
real_t qpOASES_getCPUtime( )
{
	real_t current_time = -1.0;

	#if defined(__WIN32__) || defined(WIN32)
	LARGE_INTEGER counter, frequency;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&counter);
	current_time = ((real_t) counter.QuadPart) / ((real_t) frequency.QuadPart);
	#elif defined(LINUX) || defined(__LINUX__)
	struct timeval theclock;
	gettimeofday( &theclock,0 );
	current_time = 1.0*theclock.tv_sec + 1.0e-6*theclock.tv_usec;
	#endif

	return current_time;
}



/*
 *	g e t N o r m
 */
real_t qpOASES_getNorm( const real_t* const v, int n, int type )
{
	int i;

	real_t norm = 0.0;

	switch ( type )
	{
		case 2:
			for( i=0; i<n; ++i )
				norm += v[i]*v[i];
			return qpOASES_getSqrt( norm );

		case 1:
			for( i=0; i<n; ++i )
				norm += qpOASES_getAbs( v[i] );
			return norm;

		default:
			THROWERROR( RET_INVALID_ARGUMENTS );
			return -QPOASES_INFTY;
	}
}



/*
 *	g e t K k t V i o l a t i o n
 */
returnValue qpOASES_getKktViolation(	int nV, int nC,
										const real_t* const H, const real_t* const g,
										const real_t* const A, const real_t* const lb, const real_t* const ub,
										const real_t* const lbA, const real_t* const ubA,
										const real_t* const x, const real_t* const y,
										real_t* const _stat, real_t* const feas, real_t* const cmpl
										)
{
	/* Tolerance for dual variables considered zero. */
	const real_t dualActiveTolerance = 1.0e3 * QPOASES_EPS;

	int i, j;
	real_t sum, prod;

	/* Initialize residuals */
	*_stat = *feas = *cmpl = 0.0;

	/* check stationarity */
	for (i = 0; i < nV; i++)
	{
		/* g term and variable bounds dual term */
		if ( g != 0 )
			sum = g[i] - y[i];
		else
			sum = 0.0 - y[i];

		/* H*x term */
		if ( H != 0 )
			for (j = 0; j < nV; j++) sum += H[i*nV+j] * x[j];

		/* A'*y term */
		if ( A != 0 )
			for (j = 0; j < nC; j++) sum -= A[j*nV+i] * y[nV+j];

		/* update stat */
		if (qpOASES_getAbs(sum) > *_stat)
			*_stat = qpOASES_getAbs(sum);
	}

	/* check primal feasibility and complementarity */
	/* variable bounds */
	for (i = 0; i < nV; i++)
	{
		/* feasibility */
		if ( lb != 0 )
			if (lb[i] - x[i] > *feas)
				*feas = lb[i] - x[i];

		if ( ub != 0 )
			if (x[i] - ub[i] > *feas)
				*feas = x[i] - ub[i];

		/* complementarity */
		prod = 0.0;

		if ( lb != 0 )
			if (y[i] > dualActiveTolerance) /* lower bound */
				prod = (x[i] - lb[i]) * y[i];

		if ( ub != 0 )
			if (y[i] < -dualActiveTolerance) /* upper bound */
				prod = (x[i] - ub[i]) * y[i];

		if (qpOASES_getAbs(prod) > *cmpl)
			*cmpl = qpOASES_getAbs(prod);
	}
	/* A*x bounds */
	for (i = 0; i < nC; i++)
	{
		/* compute sum = (A*x)_i */
		sum = 0.0;

		if ( A != 0 )
			for (j = 0; j < nV; j++)
				sum += A[i*nV+j] * x[j];

		/* feasibility */
		if ( lbA != 0 )
			if (lbA[i] - sum > *feas)
				*feas = lbA[i] - sum;
		if ( ubA != 0 )
			if (sum - ubA[i] > *feas)
				*feas = sum - ubA[i];

		/* complementarity */
		prod = 0.0;

		if ( lbA != 0 )
			if (y[nV+i] > dualActiveTolerance) /* lower bound */
				prod = (sum - lbA[i]) * y[nV+i];

		if ( ubA != 0 )
			if (y[nV+i] < -dualActiveTolerance) /* upper bound */
				prod = (sum - ubA[i]) * y[nV+i];

		if (qpOASES_getAbs(prod) > *cmpl)
			*cmpl = qpOASES_getAbs(prod);
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t K k t V i o l a t i o n S B
 */
returnValue qpOASES_getKktViolationSB(	int nV,
										const real_t* const H, const real_t* const g,
										const real_t* const lb, const real_t* const ub,
										const real_t* const x, const real_t* const y,
										real_t* const _stat, real_t* const feas, real_t* const cmpl
										)
{
	return qpOASES_getKktViolation(	nV,0,
									H,g,0,lb,ub,0,0,
									x,y,
									_stat,feas,cmpl
									);
}


/*
 *	c o n v e r t B o o l e a n T y p e T o S t r i n g
 */
returnValue qpOASES_convertBooleanTypeToString( BooleanType value, char* const string )
{
	#ifndef __SUPPRESSANYOUTPUT__
	if ( value == BT_FALSE )
		snprintf( string,20,"BT_FALSE" );
	else
		snprintf( string,20,"BT_TRUE" );
	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	c o n v e r t S u b j e c t T o S t a t u s T o S t r i n g
 */
returnValue qpOASES_convertSubjectToStatusToString( SubjectToStatus value, char* const string )
{
	#ifndef __SUPPRESSANYOUTPUT__
	switch( value )
	{
		case ST_INACTIVE:
			snprintf( string,20,"ST_INACTIVE" );
			break;

		case ST_LOWER:
			snprintf( string,20,"ST_LOWER" );
			break;

		case ST_UPPER:
			snprintf( string,20,"ST_UPPER" );
			break;

		case ST_UNDEFINED:
			snprintf( string,20,"ST_UNDEFINED" );
			break;

		case ST_INFEASIBLE_LOWER:
			snprintf( string,20,"ST_INFEASIBLE_LOWER" );
			break;

		case ST_INFEASIBLE_UPPER:
			snprintf( string,20,"ST_INFEASIBLE_UPPER" );
			break;

		default:
			snprintf( string,20,"<invalid value>" );
			break;
	}
	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	c o n v e r t P r i n t L  e v e l T o S t r i n g
 */
returnValue qpOASES_convertPrintLevelToString( PrintLevel value, char* const string )
{
	#ifndef __SUPPRESSANYOUTPUT__
	switch( value )
	{
		case PL_NONE:
			snprintf( string,20,"PL_NONE" );
			break;

		case PL_LOW:
			snprintf( string,20,"PL_LOW" );
			break;

		case PL_MEDIUM:
			snprintf( string,20,"PL_MEDIUM" );
			break;

		case PL_HIGH:
			snprintf( string,20,"PL_HIGH" );
			break;

		case PL_TABULAR:
			snprintf( string,20,"PL_TABULAR" );
			break;

		case PL_DEBUG_ITER:
			snprintf( string,20,"PL_DEBUG_ITER" );
			break;

		default:
			snprintf( string,20,"<invalid value>" );
			break;
	}
	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t S i m p l e S t a t u s
 */
int qpOASES_getSimpleStatus(	returnValue returnvalue,
								BooleanType doPrintStatus
								)
{
	int retValNumber;
	int simpleStatus = -1;

	/* determine simple status from returnvalue */
	switch ( returnvalue )
	{
		case SUCCESSFUL_RETURN:
			simpleStatus = 0;
			break;

		case RET_MAX_NWSR_REACHED:
			simpleStatus = 1;
			break;

		case RET_INIT_FAILED_INFEASIBILITY:
		case RET_HOTSTART_STOPPED_INFEASIBILITY:
			simpleStatus = -2;
			break;

		case RET_INIT_FAILED_UNBOUNDEDNESS:
		case RET_HOTSTART_STOPPED_UNBOUNDEDNESS:
			simpleStatus = -3;
			break;

		default:
			simpleStatus = -1;
			break;
	}

	if ( doPrintStatus == BT_TRUE )
	{
		VisibilityStatus vsInfo = MessageHandling_getInfoVisibilityStatus( qpOASES_getGlobalMessageHandler() );
		MessageHandling_setInfoVisibilityStatus( qpOASES_getGlobalMessageHandler(),VS_VISIBLE );
		MessageHandling_setErrorCount( qpOASES_getGlobalMessageHandler(),-1 );

		retValNumber = (int)RET_SIMPLE_STATUS_P0 - simpleStatus;
		THROWINFO( (returnValue)retValNumber );

		MessageHandling_setInfoVisibilityStatus( qpOASES_getGlobalMessageHandler(),vsInfo );
	}

	return simpleStatus;
}


/*
 *	n o r m a l i s e C o n s t r a i n t s
 */
returnValue normaliseConstraints(	int nV, int nC,
									real_t* A, real_t* lbA, real_t* ubA,
									int type
									)
{
	int ii, jj;
	real_t curNorm;

	if ( ( nV <= 0 ) || ( nC <= 0 ) || ( A == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	for( ii=0; ii<nC; ++ii )
	{
		/* get row norm */
		curNorm = qpOASES_getNorm( &(A[ii*nV]),nV,type );

		if ( curNorm > QPOASES_EPS )
		{
			/* normalise if norm is positive */
			for( jj=0; jj<nV; ++jj )
				A[ii*nV + jj] /= curNorm;

			if ( lbA != 0 ) lbA[ii] /= curNorm;
			if ( ubA != 0 ) ubA[ii] /= curNorm;
		}
		else
		{
			/* if row norm is (close to) zero, kind of erase constraint */
			if ( type == 1 )
			{
				for( jj=0; jj<nV; ++jj )
					A[ii*nV + jj] = 1.0 / ((real_t)nV);
			}
			else
			{
				/* assume type == 2 */
				for( jj=0; jj<nV; ++jj )
					A[ii*nV + jj] = 1.0 / qpOASES_getSqrt((real_t)nV);
			}

			if ( lbA != 0 ) lbA[ii] = -QPOASES_INFTY;
			if ( ubA != 0 ) ubA[ii] =  QPOASES_INFTY;
		}
	}

	return SUCCESSFUL_RETURN;
}


#ifdef __DEBUG__
/*
 *	g d b _ p r i n t m at
 */
void gdb_printmat( const char *fname, real_t *M, int n, int m, int ldim )
{
	#ifdef __SUPPRESSANYOUTPUT__

	int i, j;
	FILE *fid;

	fid = fopen(fname, "wt");
	if (!fid)
	{
		perror("Error opening file: ");
		return;
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
			fprintf(fid, " %23.16e", M[j*ldim+i]);
		fprintf(fid, "\n");
	}
	fclose(fid);

	#endif /* __SUPPRESSANYOUTPUT__ */
}
#endif /* __DEBUG__ */


#if defined(__DSPACE__) || defined(__XPCTARGET__)
/*
 *	_ _ c x a _ p u r e _ v i r t u a l
 */
void __cxa_pure_virtual( void )
{
	/* put your customized implementation here! */
}
#endif /* __DSPACE__ || __XPCTARGET__*/


END_NAMESPACE_QPOASES


/*
 *	end of file
 */

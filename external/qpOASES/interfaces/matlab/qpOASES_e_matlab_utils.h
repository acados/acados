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
 *	\file interfaces/matlab/qpOASES_e_matlab_utils.h
 *	\author Hans Joachim Ferreau, Alexander Buchner
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Collects utility functions for Interface to Matlab(R) that
 *	enables to call qpOASES as a MEX function.
 *
 */



/* Work-around for settings where mexErrMsgTxt causes unexpected behaviour. */
#ifdef __AVOID_MEXERRMSGTXT__
	#define myMexErrMsgTxt( TEXT ) mexPrintf( "%s\n\n",(TEXT) );
#else
	#define myMexErrMsgTxt mexErrMsgTxt
#endif


#include "mex.h"
#include "matrix.h"
#include "string.h"



static int QPInstance_nexthandle = -1;

/*
 * QProblem instance class
 */
typedef struct
{
	int handle;

	QProblem  sqp;
	QProblemB qpb;
	BooleanType isSimplyBounded;

	DenseMatrix H;
	DenseMatrix A;

} QPInstance;


void QPInstanceCON(	QPInstance* _THIS,
					int _nV,
					int _nC,
					HessianType _hessianType,
					BooleanType _isSimplyBounded
					);
	
int QPInstance_getNV( QPInstance* _THIS );

int QPInstance_getNC( QPInstance* _THIS );


#define MAX_NUM_QPINSTANCES 10
static QPInstance QPInstances[MAX_NUM_QPINSTANCES];


/*
 *	end of file
 */

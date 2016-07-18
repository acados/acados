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
 *	\file interfaces/simulink/qpOASES_e_QProblem.c
 *	\author Hans Joachim Ferreau (thanks to Aude Perrin)
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Interface for Simulink(R) that enables to call qpOASES as a S function
 *  (variant for QPs with fixed matrices).
 *
 */


#include <stdlib.h>

#include <qpOASES_e.h>
#include "qpOASES_e_simulink_utils.c"


#ifdef __cplusplus
extern "C" {
#endif


#define S_FUNCTION_NAME   qpOASES_e_QProblem	/**< Name of the S function. */
#define S_FUNCTION_LEVEL  2						/**< S function level. */

#define MDL_START								/**< Activate call to mdlStart. */

#include "simstruc.h"


/* SETTINGS: */
#define SAMPLINGTIME    -1						/**< Sampling time. */
#define NCONTROLINPUTS  2						/**< Number of control inputs. */
#define MAXITER         100						/**< Maximum number of iterations. */
#define HESSIANTYPE     HST_UNKNOWN				/**< Hessian type, see documentation of QProblem class constructor. */



static void mdlInitializeSizes (SimStruct *S)   /* Init sizes array */
{
	int nU = NCONTROLINPUTS;

	/* Specify the number of continuous and discrete states */
	ssSetNumContStates(S, 0);
	ssSetNumDiscStates(S, 0);

	/* Specify the number of parameters */
	ssSetNumSFcnParams(S, 2); /* H, A */
	if ( ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S) )
		return;

	/* Specify the number of intput ports */
	if ( !ssSetNumInputPorts(S, 5) )
		return;

	/* Specify the number of output ports */
	if ( !ssSetNumOutputPorts(S, 4) )
		return;

	/* Specify dimension information for the input ports */
	ssSetInputPortVectorDimension(S, 0, DYNAMICALLY_SIZED); /* g */
	ssSetInputPortVectorDimension(S, 1, DYNAMICALLY_SIZED); /* lb */
	ssSetInputPortVectorDimension(S, 2, DYNAMICALLY_SIZED); /* ub */
	ssSetInputPortVectorDimension(S, 3, DYNAMICALLY_SIZED); /* lbA */
	ssSetInputPortVectorDimension(S, 4, DYNAMICALLY_SIZED); /* ubA */

	/* Specify dimension information for the output ports */
	ssSetOutputPortVectorDimension(S, 0, nU );  /* uOpt */
	ssSetOutputPortVectorDimension(S, 1, 1 );   /* fval */
	ssSetOutputPortVectorDimension(S, 2, 1 );   /* exitflag */
	ssSetOutputPortVectorDimension(S, 3, 1 );   /* iter */

	/* Specify the direct feedthrough status */
	ssSetInputPortDirectFeedThrough(S, 0, 1);
	ssSetInputPortDirectFeedThrough(S, 1, 1);
	ssSetInputPortDirectFeedThrough(S, 2, 1);
	ssSetInputPortDirectFeedThrough(S, 3, 1);
	ssSetInputPortDirectFeedThrough(S, 4, 1);
	ssSetInputPortDirectFeedThrough(S, 5, 1);
	ssSetInputPortDirectFeedThrough(S, 6, 1);

	/* One sample time */
	ssSetNumSampleTimes(S, 1);

	/* global variables:
     * 0: problem
     * 1: H
     * 2: g
     * 3: A
     * 4: lb
     * 5: ub
     * 6: lbA
     * 7: ubA
     */

	/* Specify the size of the block's pointer work vector */
    ssSetNumPWork(S, 8);
}


#if defined(MATLAB_MEX_FILE)

#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
	if ( !ssSetInputPortDimensionInfo(S, port, dimsInfo) )
		return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
	if ( !ssSetOutputPortDimensionInfo(S, port, dimsInfo) )
		return;
}

#endif


static void mdlInitializeSampleTimes(SimStruct *S)
{
	ssSetSampleTime(S, 0, SAMPLINGTIME);
	ssSetOffsetTime(S, 0, 0.0);
}


static void mdlStart(SimStruct *S)
{
	USING_NAMESPACE_QPOASES

	int nU = NCONTROLINPUTS;
	int size_g, size_lb, size_ub, size_lbA, size_ubA;
	int size_H, nRows_H, nCols_H, size_A, nRows_A, nCols_A;
	int nV, nC;

	static QProblem problem;
	static Options problemOptions;


	/* get block inputs dimensions */
	const mxArray* in_H = ssGetSFcnParam(S, 0);
	const mxArray* in_A = ssGetSFcnParam(S, 1);

	if ( mxIsEmpty(in_H) == 1 )
	{
		if ( ( HESSIANTYPE != HST_ZERO ) && ( HESSIANTYPE != HST_IDENTITY ) )
		{
			#ifndef __SUPPRESSANYOUTPUT__
			mexErrMsgTxt( "ERROR (qpOASES): Hessian can only be empty if type is set to HST_ZERO or HST_IDENTITY!" );
			#endif
			return;
		}
		
	    nRows_H = 0;
		nCols_H = 0;
		size_H  = 0;
	}
	else
	{
	    nRows_H = (int)mxGetM(in_H);
		nCols_H = (int)mxGetN(in_H);
		size_H  = nRows_H * nCols_H;
	}

	if ( mxIsEmpty(in_A) == 1 )
	{
	    nRows_A = 0;
		nCols_A = 0;
		size_A  = 0;
	}
	else
	{
	    nRows_A = (int)mxGetM(in_A);
		nCols_A = (int)mxGetN(in_A);
		size_A  = nRows_A * nCols_A;
	}

	size_g   = ssGetInputPortWidth(S, 0);
	size_lb  = ssGetInputPortWidth(S, 1);
	size_ub  = ssGetInputPortWidth(S, 2);
	size_lbA = ssGetInputPortWidth(S, 3);
	size_ubA = ssGetInputPortWidth(S, 4);


	/* dimension checks */
	nV = size_g;
	nC = nRows_A;


	if ( MAXITER < 0 )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Maximum number of iterations must not be negative!" );
		#endif
		return;
	}

	if ( nV <= 0 )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Dimension mismatch (nV must be positive)!" );
		#endif
		return;
	}

	if ( ( size_H != nV*nV ) && ( size_H != 0 ) )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Dimension mismatch in H!" );
		#endif
		return;
	}

	if ( nRows_H != nCols_H )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Hessian matrix must be square matrix!" );
		#endif
		return;
	}

	if ( ( nU < 1 ) || ( nU > nV ) )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Invalid number of control inputs!" );
		#endif
		return;
	}

	if ( ( size_lb != nV ) && ( size_lb != 0 ) )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Dimension mismatch in lb!" );
		#endif
		return;
	}

	if ( ( size_ub != nV ) && ( size_lb != 0 ) )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Dimension mismatch in ub!" );
		#endif
		return;
	}

	if ( ( size_lbA != nC ) && ( size_lbA != 0 ) )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Dimension mismatch in lbA!" );
		#endif
		return;
	}

	if ( ( size_ubA != nC ) && ( size_ubA != 0 ) )
	{
		#ifndef __SUPPRESSANYOUTPUT__
		mexErrMsgTxt( "ERROR (qpOASES): Dimension mismatch in ubA!" );
		#endif
		return;
	}


	/* allocate QProblem object */
	QProblemCON( &problem,nV,nC,HESSIANTYPE );

	Options_setToMPC( &problemOptions );
	QProblem_setOptions( &problem,problemOptions );
	#ifndef __DEBUG__
	QProblem_setPrintLevel( &problem,PL_LOW );
	#endif
	#ifdef __SUPPRESSANYOUTPUT__
	QProblem_setPrintLevel( &problem,PL_NONE );
	#endif

	ssGetPWork(S)[0] = (void*)(&problem);

	/* allocate memory for QP data ... */
	if ( size_H > 0 )
		ssGetPWork(S)[1] = (void *) calloc( size_H, sizeof(real_t) );	/* H */
	else
		ssGetPWork(S)[1] = 0;

	ssGetPWork(S)[2] = (void *) calloc( size_g, sizeof(real_t) );		/* g */
	ssGetPWork(S)[3] = (void *) calloc( size_A, sizeof(real_t) );		/* A */

	if ( size_lb > 0 )
		ssGetPWork(S)[4] = (void *) calloc( size_lb, sizeof(real_t) );	/* lb */
	else
		ssGetPWork(S)[4] = 0;

	if ( size_ub > 0 )
		ssGetPWork(S)[5] = (void *) calloc( size_ub, sizeof(real_t) );	/* ub */
	else
		ssGetPWork(S)[5] = 0;
	
	if ( size_lbA > 0 )
		ssGetPWork(S)[6] = (void *) calloc( size_lbA, sizeof(real_t) );	/* lbA */
	else
		ssGetPWork(S)[6] = 0;

	if ( size_ubA > 0 )
		ssGetPWork(S)[7] = (void *) calloc( size_ubA, sizeof(real_t) );	/* ubA */
	else
		ssGetPWork(S)[7] = 0;
}


static void mdlOutputs(SimStruct *S, int_T tid)
{
	USING_NAMESPACE_QPOASES

	int i;
	int nV, nC;
	returnValue status;

	int nWSR = MAXITER;
	int nU   = NCONTROLINPUTS;

	InputRealPtrsType in_g, in_lb, in_ub, in_lbA, in_ubA;

	QProblem* problem;
	real_t *H, *g, *A, *lb, *ub, *lbA, *ubA;

	real_t xOpt[NVMAX];
    real_t yOpt[NVMAX+NCMAX];
    
	real_T *out_uOpt, *out_objVal, *out_status, *out_nWSR;
    real_t stat, feas, cmpl, kktTol;

	int nWSR_retry;
    
    FILE* matFile = 0;
    char fileName[20] = "qpData000000000.mat";
    

	/* get pointers to block inputs ... */
	const mxArray* in_H = ssGetSFcnParam(S, 0);
	const mxArray* in_A = ssGetSFcnParam(S, 1);
	in_g   = ssGetInputPortRealSignalPtrs(S, 0);
	in_lb  = ssGetInputPortRealSignalPtrs(S, 1);
	in_ub  = ssGetInputPortRealSignalPtrs(S, 2);
	in_lbA = ssGetInputPortRealSignalPtrs(S, 3);
	in_ubA = ssGetInputPortRealSignalPtrs(S, 4);


	/* ... and to the QP data */
	problem = (QProblem*)(ssGetPWork(S)[0]);

	H   = (real_t*)(ssGetPWork(S)[1]);
	g   = (real_t*)(ssGetPWork(S)[2]);
	A   = (real_t*)(ssGetPWork(S)[3]);
	lb  = (real_t*)(ssGetPWork(S)[4]);
	ub  = (real_t*)(ssGetPWork(S)[5]);
	lbA = (real_t*)(ssGetPWork(S)[6]);
	ubA = (real_t*)(ssGetPWork(S)[7]);


	/* setup QP data */
	nV = ssGetInputPortWidth(S, 0); /* nV = size_g */
	nC = (int)mxGetM(in_A);			/* nC = nRows_A*/

	if ( H != 0 )
	{
		/* no conversion from FORTRAN to C as Hessian is symmetric! */
		for ( i=0; i<nV*nV; ++i )
			H[i] = (mxGetPr(in_H))[i];
	}

	convertFortranToC( mxGetPr(in_A),nV,nC, A );

	for ( i=0; i<nV; ++i )
		g[i] = (*in_g)[i];

	if ( lb != 0 )
	{
		for ( i=0; i<nV; ++i )
			lb[i] = (*in_lb)[i];
	}

	if ( ub != 0 )
	{
		for ( i=0; i<nV; ++i )
			ub[i] = (*in_ub)[i];
	}

	if ( lbA != 0 )
	{
		for ( i=0; i<nC; ++i )
			lbA[i] = (*in_lbA)[i];
	}

	if ( ubA != 0 )
	{
		for ( i=0; i<nC; ++i )
			ubA[i] = (*in_ubA)[i];
	}
    
   
    #ifdef __SIMULINK_DEBUG__

    sprintf( fileName,"qpData%09d.mat",(int)(QProblem_getCount( problem )) );
    matFile = fopen( fileName,"w+" );
    /*qpOASES_writeIntoMatFile( matFile, H,   nV,nV, "H"   );*/
    //qpOASES_writeIntoMatFile( matFile, g,   nV,1,  "g"   );
    qpOASES_writeIntoMatFile( matFile, lb,  nV,1,  "lb"  );
    qpOASES_writeIntoMatFile( matFile, ub,  nV,1,  "ub"  );
    //qpOASES_writeIntoMatFile( matFile, A,   nC,nV, "A"   );
    qpOASES_writeIntoMatFile( matFile, lbA, nC,1,  "lbA" );
    qpOASES_writeIntoMatFile( matFile, ubA, nC,1,  "ubA" );
    fclose( matFile );

    #endif /* __SIMULINK_DEBUG__ */


	if ( QProblem_getCount( problem ) == 0 )
	{
		/* initialise and solve first QP */
		status = QProblem_init( problem,H,g,A,lb,ub,lbA,ubA, &nWSR,0 );
		QProblem_getPrimalSolution( problem,xOpt );
        QProblem_getDualSolution(   problem,yOpt );
	}
	else
	{
		/* solve neighbouring QP using hotstart technique */
		status = QProblem_hotstart( problem,g,lb,ub,lbA,ubA, &nWSR,0 );
		if ( ( status != SUCCESSFUL_RETURN ) && ( status != RET_MAX_NWSR_REACHED ) )
		{
			/* if an error occurs, reset problem data structures ... */
			QProblem_reset( problem );
            
            /* ... and initialise/solve again with remaining number of iterations. */
            nWSR_retry = MAXITER - nWSR;
			status = QProblem_init( problem,H,g,A,lb,ub,lbA,ubA, &nWSR_retry,0 );
			nWSR += nWSR_retry;
			
		}

		/* obtain optimal solution */
		QProblem_getPrimalSolution( problem,xOpt );
        QProblem_getDualSolution(   problem,yOpt );
	}

	/* generate block output: status information ... */
	out_uOpt   = ssGetOutputPortRealSignal(S, 0);
	out_objVal = ssGetOutputPortRealSignal(S, 1);
	out_status = ssGetOutputPortRealSignal(S, 2);
	out_nWSR   = ssGetOutputPortRealSignal(S, 3);

	for ( i=0; i<nU; ++i )
		out_uOpt[i] = (real_T)(xOpt[i]);

	out_objVal[0] = (real_T)(QProblem_getObjVal( problem ));
	out_status[0] = (real_t)(qpOASES_getSimpleStatus( status,BT_FALSE ));
	out_nWSR[0]   = (real_T)(nWSR);
    
	removeNaNs( out_uOpt,nU );
	removeInfs( out_uOpt,nU );
	removeNaNs( out_objVal,1 );
	removeInfs( out_objVal,1 );

    #ifdef __SIMULINK_DEBUG__

    qpOASES_getKktViolation(	nV,nC,
								H,g,A,lb,ub,lbA,ubA,
								xOpt,yOpt,
								&stat,&feas,&cmpl
								);
    mexPrintf( "KKT residuals:  stat=%e, feas=%e, cmpl=%e\n", stat,feas,cmpl );

	kktTol = stat;
	if ( feas > kktTol )
		kktTol = feas;
	if ( cmpl > kktTol )
		kktTol = cmpl;

	if ( ( qpOASES_getAbs( out_status[0] ) < QPOASES_EPS ) && ( kktTol > 1e-4 ) )
		out_status[0] = 2.0;

    #endif /* __SIMULINK_DEBUG__ */
}


static void mdlTerminate(SimStruct *S)
{
	USING_NAMESPACE_QPOASES

	int i;

	/* reset global message handler */
	MessageHandling_reset( qpOASES_getGlobalMessageHandler() );

	for ( i=1; i<8; ++i )
	{
		if ( ssGetPWork(S)[i] != 0 )
			free( ssGetPWork(S)[i] );
	}
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif


#ifdef __cplusplus
}
#endif


/*
 *	end of file
 */

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
 *	\file src/OQPinterface.cpp
 *	\author Hans Joachim Ferreau
 *	\version 3.1embedded
 *	\date 2008-2015
 *
 *	Implementation of an interface comprising several utility functions
 *	for solving test problems from the Online QP Benchmark Collection
 *	(This collection is no longer maintained, see 
 *	http://www.qpOASES.org/onlineQP for a backup).
 *
 */


#include <qpOASES_e/extras/OQPinterface.h>
#include <qpOASES_e/QProblemB.h>
#include <qpOASES_e/QProblem.h>


BEGIN_NAMESPACE_QPOASES


/*
 *	r e a d O Q P d i m e n s i o n s
 */
returnValue readOQPdimensions(	const char* path,
								int* nQP, int* nV, int* nC, int* nEC
								)
{
	int dims[4];
	
	/* 1) Setup file name where dimensions are stored. */
	char filename[QPOASES_MAX_STRING_LENGTH];
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sdims.oqp",path );

	/* 2) Load dimensions from file. */
	if ( qpOASES_readFromFileI( dims,4,filename ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	*nQP = dims[0];
	*nV  = dims[1];
	*nC  = dims[2];
	*nEC = dims[3];

	/* printf( "nQP = %d,  nV = %d,  nC = %d,  nEC = %d\n",*nQP,*nV,*nC,*nEC ); */

	/* consistency check */
	if ( ( *nQP <= 0 ) || ( *nV <= 0 ) || ( *nC < 0 ) || ( *nEC < 0 ) )
		return THROWERROR( RET_FILEDATA_INCONSISTENT );

	if ( ( *nV > NVMAX ) || ( *nC > NCMAX ) || ( *nQP > NQPMAX ) )
		return THROWERROR( RET_UNABLE_TO_READ_BENCHMARK );

	return SUCCESSFUL_RETURN;
}


/*
 *	r e a d O Q P d a t a
 */
returnValue readOQPdata(	const char* path,
							int* nQP, int* nV, int* nC, int* nEC,
							real_t* H, real_t* g, real_t* A, real_t* lb, real_t* ub, real_t* lbA, real_t* ubA,
							real_t* xOpt, real_t* yOpt, real_t* objOpt
							)
{
	char filename[QPOASES_MAX_STRING_LENGTH];

	/* consistency check */
	if ( ( H == 0 ) || ( g == 0 ) || ( lb == 0 ) || ( ub == 0 ) )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Obtain OQP dimensions. */
	if ( readOQPdimensions( path, nQP,nV,nC,nEC ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );


	/* another consistency check */
	if ( ( *nC > 0 ) && ( ( A == 0 ) || ( lbA == 0 ) || ( ubA == 0 ) ) )
		return THROWERROR( RET_FILEDATA_INCONSISTENT );


	/* 2) Allocate memory and load OQP data: */
	/* Hessian matrix */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sH.oqp",path );
	if ( qpOASES_readFromFileM( H,(*nV),(*nV),filename ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* gradient vector sequence */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sg.oqp",path );
	if ( qpOASES_readFromFileM( g,(*nQP),(*nV),filename ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	
	/* lower bound vector sequence */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%slb.oqp",path );
	if ( qpOASES_readFromFileM( lb,(*nQP),(*nV),filename ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );

	/* upper bound vector sequence */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sub.oqp",path );
	if ( qpOASES_readFromFileM( ub,(*nQP),(*nV),filename ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_UNABLE_TO_READ_FILE );
	
	if ( (*nC) > 0 )
	{
		/* Constraint matrix */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sA.oqp",path );
		if ( qpOASES_readFromFileM( A,(*nC),(*nV),filename ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_UNABLE_TO_READ_FILE );

		/* lower constraints' bound vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%slbA.oqp",path );
		if ( qpOASES_readFromFileM( lbA,(*nQP),(*nC),filename ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_UNABLE_TO_READ_FILE );

		/* upper constraints' bound vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%subA.oqp",path );
		if ( qpOASES_readFromFileM( ubA,(*nQP),(*nC),filename ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	if ( xOpt != 0 )
	{
		/* primal solution vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sx_opt.oqp",path );
		if ( qpOASES_readFromFileM( xOpt,(*nQP),(*nV),filename ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	if ( yOpt != 0 )
	{
		/* dual solution vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sy_opt.oqp",path );
		if ( qpOASES_readFromFileM( yOpt,(*nQP),(*nV)+(*nC),filename ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}

	if ( objOpt != 0 )
	{
		/* dual solution vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sobj_opt.oqp",path );
		if ( qpOASES_readFromFileM( objOpt,(*nQP),1,filename ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_UNABLE_TO_READ_FILE );
	}
	
	return SUCCESSFUL_RETURN;
}


/*
 *	s o l v e O Q P b e n c h m a r k
 */
returnValue solveOQPbenchmark(	int nQP, int nV, int nC, int nEC,
								real_t* _H, const real_t* const g, real_t* _A,
								const real_t* const lb, const real_t* const ub,
								const real_t* const lbA, const real_t* const ubA,
								BooleanType isSparse, BooleanType useHotstarts, 
								const Options* options, int maxAllowedNWSR,
								real_t* maxNWSR, real_t* avgNWSR, real_t* maxCPUtime, real_t* avgCPUtime,
								real_t* maxStationarity, real_t* maxFeasibility, real_t* maxComplementarity
								)
{
	int k;
	
	myStatic QProblem qp;
	returnValue returnvalue;

	/* I) SETUP AUXILIARY VARIABLES: */
	/* 1) Keep nWSR and store current and maximum number of
	 *    working set recalculations in temporary variables */
	int nWSRcur;

	real_t CPUtimeLimit = *maxCPUtime;
	real_t CPUtimeCur = CPUtimeLimit;
	real_t stat, feas, cmpl;
	
	/* 2) Pointers to data of current QP ... */
	const real_t* gCur;
	const real_t* lbCur;
	const real_t* ubCur;
	const real_t* lbACur;
	const real_t* ubACur;

	/* 3) Vectors for solution obtained by qpOASES. */
	myStatic real_t x[NVMAX];
	myStatic real_t y[NVMAX+NCMAX];

	/* 4) Prepare matrix objects */
	DenseMatrix *H, *A;
	myStatic DenseMatrix HH, AA;


	DenseMatrixCON( &HH, nV, nV, nV, _H );
	DenseMatrixCON( &AA, nC, nV, nV, _A );
	
	H = &HH;
	A = &AA;
	
	*maxNWSR = 0;
	*avgNWSR = 0;
	*maxCPUtime = 0.0;
	*avgCPUtime = 0.0;
	*maxStationarity = 0.0;
	*maxFeasibility = 0.0;
	*maxComplementarity = 0.0;
	
	/*DenseMatrix_print( H );*/

	/* II) SETUP QPROBLEM OBJECT */
	QProblemCON( &qp,nV,nC,HST_UNKNOWN );
	QProblem_setOptions( &qp,*options );
	/*QProblem_setPrintLevel( &qp,PL_LOW );*/

	 /* QProblem_printOptions( &qp ); */

	/* III) RUN BENCHMARK SEQUENCE: */

	for( k=0; k<nQP; ++k )
	{
		/* 1) Update pointers to current QP data. */
		gCur   = &( g[k*nV] );
		lbCur  = &( lb[k*nV] );
		ubCur  = &( ub[k*nV] );
		lbACur = &( lbA[k*nC] );
		ubACur = &( ubA[k*nC] );

		/* 2) Set nWSR and maximum CPU time. */
		nWSRcur = maxAllowedNWSR;
		CPUtimeCur = CPUtimeLimit;

		/* 3) Solve current QP. */
		if ( ( k == 0 ) || ( useHotstarts == BT_FALSE ) )
		{
			/* initialise */
			returnvalue = QProblem_initM( &qp, H,gCur,A,lbCur,ubCur,lbACur,ubACur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* hotstart */
			returnvalue = QProblem_hotstart( &qp, gCur,lbCur,ubCur,lbACur,ubACur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}

		/* 4) Obtain solution vectors and objective function value */
		QProblem_getPrimalSolution( &qp,x );
		QProblem_getDualSolution( &qp,y );

		/* 5) Compute KKT residuals */
		qpOASES_getKktViolation( nV,nC, _H,gCur,_A,lbCur,ubCur,lbACur,ubACur, x,y, &stat,&feas,&cmpl );
		
		/* 6) Update maximum values. */
		if ( nWSRcur > *maxNWSR )
			*maxNWSR = nWSRcur;
		if (stat > *maxStationarity) *maxStationarity = stat;
		if (feas > *maxFeasibility) *maxFeasibility = feas;
		if (cmpl > *maxComplementarity) *maxComplementarity = cmpl;

		if ( CPUtimeCur > *maxCPUtime )
			*maxCPUtime = CPUtimeCur;
	
		*avgNWSR += nWSRcur;
		*avgCPUtime += CPUtimeCur;
	}
	*avgNWSR /= nQP;
	*avgCPUtime /= ((double)nQP);

	return SUCCESSFUL_RETURN;
}


/*
 *	s o l v e O Q P b e n c h m a r k
 */
returnValue solveOQPbenchmarkB(	int nQP, int nV,
								real_t* _H, const real_t* const g,
								const real_t* const lb, const real_t* const ub,
								BooleanType isSparse, BooleanType useHotstarts, 
								const Options* options, int maxAllowedNWSR,
								real_t* maxNWSR, real_t* avgNWSR, real_t* maxCPUtime, real_t* avgCPUtime,
								real_t* maxStationarity, real_t* maxFeasibility, real_t* maxComplementarity
								)
{
	int k;
	
	myStatic QProblemB qp;
	returnValue returnvalue;

	/* I) SETUP AUXILIARY VARIABLES: */
	/* 1) Keep nWSR and store current and maximum number of
	 *    working set recalculations in temporary variables */
	int nWSRcur;

	real_t CPUtimeLimit = *maxCPUtime;
	real_t CPUtimeCur = CPUtimeLimit;
	real_t stat, feas, cmpl;

	/* 2) Pointers to data of current QP ... */
	const real_t* gCur;
	const real_t* lbCur;
	const real_t* ubCur;

	/* 3) Vectors for solution obtained by qpOASES. */
	myStatic real_t x[NVMAX];
	myStatic real_t y[NVMAX];

	/* 4) Prepare matrix objects */
	DenseMatrix *H;
	myStatic DenseMatrix HH;
	
	DenseMatrixCON( &HH, nV, nV, nV, _H );
	H = &HH;
	
	*maxNWSR = 0;
	*avgNWSR = 0;
	*maxCPUtime = 0.0;
	*avgCPUtime = 0.0;
	*maxStationarity = 0.0;
	*maxFeasibility = 0.0;
	*maxComplementarity = 0.0;

	/* II) SETUP QPROBLEM OBJECT */
	QProblemBCON( &qp,nV,HST_UNKNOWN );
	QProblemB_setOptions( &qp,*options );
	/*QProblemB_setPrintLevel( &qp,PL_LOW );*/


	/* III) RUN BENCHMARK SEQUENCE: */
	for( k=0; k<nQP; ++k )
	{
		/* 1) Update pointers to current QP data. */
		gCur   = &( g[k*nV] );
		lbCur  = &( lb[k*nV] );
		ubCur  = &( ub[k*nV] );

		/* 2) Set nWSR and maximum CPU time. */
		nWSRcur = maxAllowedNWSR;
		CPUtimeCur = CPUtimeLimit;

		/* 3) Solve current QP. */
		if ( ( k == 0 ) || ( useHotstarts == BT_FALSE ) )
		{
			/* initialise */
			returnvalue = QProblemB_initM( &qp,H,gCur,lbCur,ubCur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* hotstart */
			returnvalue = QProblemB_hotstart( &qp,gCur,lbCur,ubCur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}

		/* 4) Obtain solution vectors and objective function value ... */
		QProblemB_getPrimalSolution( &qp,x );
		QProblemB_getDualSolution( &qp,y );

		/* 5) Compute KKT residuals */
		qpOASES_getKktViolationSB( nV, _H,gCur,lbCur,ubCur, x,y, &stat,&feas,&cmpl );

		/* 6) update maximum values. */
		if ( nWSRcur > *maxNWSR )
			*maxNWSR = nWSRcur;
		if (stat > *maxStationarity) *maxStationarity = stat;
		if (feas > *maxFeasibility) *maxFeasibility = feas;
		if (cmpl > *maxComplementarity) *maxComplementarity = cmpl;

		if ( CPUtimeCur > *maxCPUtime )
			*maxCPUtime = CPUtimeCur;
		
		*avgNWSR += nWSRcur;
		*avgCPUtime += CPUtimeCur;
	}
	*avgNWSR /= nQP;
	*avgCPUtime /= ((double)nQP);

	return SUCCESSFUL_RETURN;
}


/*
 *	r u n O Q P b e n c h m a r k
 */
returnValue runOQPbenchmark(	const char* path, BooleanType isSparse, BooleanType useHotstarts, 
								const Options* options, int maxAllowedNWSR,
								real_t* maxNWSR, real_t* avgNWSR, real_t* maxCPUtime, real_t* avgCPUtime,
								real_t* maxStationarity, real_t* maxFeasibility, real_t* maxComplementarity
								)
{
	int nQP=0, nV=0, nC=0, nEC=0;

	myStatic real_t H[NVMAX*NVMAX];
	myStatic real_t g[NQPMAX*NVMAX];
	myStatic real_t A[NCMAX*NVMAX];
	myStatic real_t lb[NQPMAX*NVMAX];
	myStatic real_t ub[NQPMAX*NVMAX];
	myStatic real_t lbA[NQPMAX*NCMAX];
	myStatic real_t ubA[NQPMAX*NCMAX];

	returnValue returnvalue;

	/* I) SETUP BENCHMARK: */
	/* 1) Obtain QP sequence dimensions. */
	/*if ( readOQPdimensions( path, &nQP,&nV,&nC,&nEC ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_BENCHMARK_ABORTED );*/

	/* 2) Read OQP benchmark data. */
	if ( readOQPdata(	path,
						&nQP,&nV,&nC,&nEC,
						H,g,A,lb,ub,lbA,ubA,
						0,0,0
						) != SUCCESSFUL_RETURN )
	{
		return THROWERROR( RET_UNABLE_TO_READ_BENCHMARK );
	}
	
	/* II) SOLVE BENCHMARK */
	if ( nC > 0 )
	{
		returnvalue = solveOQPbenchmark(	nQP,nV,nC,nEC,
											H,g,A,lb,ub,lbA,ubA,
											isSparse,useHotstarts,
											options,maxAllowedNWSR,
											maxNWSR,avgNWSR,maxCPUtime,avgCPUtime,
											maxStationarity,maxFeasibility,maxComplementarity
											);

		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}
	else
	{
		returnvalue = solveOQPbenchmarkB(	nQP,nV,
											H,g,lb,ub,
											isSparse,useHotstarts,
											options,maxAllowedNWSR,
											maxNWSR,avgNWSR,maxCPUtime,avgCPUtime,
											maxStationarity,maxFeasibility,maxComplementarity
											);

		if ( returnvalue != SUCCESSFUL_RETURN )
			return THROWERROR( returnvalue );
	}

	return SUCCESSFUL_RETURN;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */

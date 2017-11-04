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


BEGIN_NAMESPACE_QPOASES

int OQPbenchmark_ws_calculateMemorySize( unsigned int nV, unsigned int nC )
{
	int size = 0;
	size += sizeof(OQPbenchmark_ws);
	size += QProblem_calculateMemorySize(nV, nC);     // qp
	size += DenseMatrix_calculateMemorySize(nV, nV);  // H
	size += DenseMatrix_calculateMemorySize(nC, nV);  // A
	size += 1 * nV * sizeof(real_t);             	  // x
	size += 1 * (nV + nC) * sizeof(real_t);			  // y

	size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
	size += 1 * 64;                // align once to typical cache line size

	return size;
}

char *OQPbenchmark_ws_assignMemory( unsigned int nV, unsigned int nC, OQPbenchmark_ws **mem, void *raw_memory )
{
	// char pointer
	char *c_ptr = (char *)raw_memory;

	// assign structures
	*mem = (OQPbenchmark_ws *) c_ptr;
	c_ptr += sizeof(OQPbenchmark_ws);

	(*mem)->qp = (QProblem *) c_ptr;
	c_ptr = QProblem_assignMemory(nV, nC, &((*mem)->qp), c_ptr);

	(*mem)->H = (DenseMatrix *) c_ptr;
	c_ptr = DenseMatrix_assignMemory(nV, nV, &((*mem)->H), c_ptr);

	(*mem)->A = (DenseMatrix *) c_ptr;
	c_ptr = DenseMatrix_assignMemory(nC, nV, &((*mem)->A), c_ptr);

	// align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
	c_ptr = (char *)s_ptr;

	// assign data
	(*mem)->x = (real_t *) c_ptr;
	c_ptr += nV * sizeof(real_t);

	(*mem)->y = (real_t *) c_ptr;
	c_ptr += (nV + nC) * sizeof(real_t);

	return c_ptr;
}

OQPbenchmark_ws *OQPbenchmark_ws_createMemory( unsigned int nV, unsigned int nC )
{
	OQPbenchmark_ws *mem;
    int memory_size = OQPbenchmark_ws_calculateMemorySize(nV, nC);
    void *raw_memory_ptr = malloc(memory_size);
    char *ptr_end =  OQPbenchmark_ws_assignMemory(nV, nC, &mem, raw_memory_ptr);
    assert((char*)raw_memory_ptr + memory_size >= ptr_end); (void) ptr_end;
    return mem;
}

int OQPbenchmarkB_ws_calculateMemorySize( unsigned int nV )
{
	int size = 0;
	size += sizeof(OQPbenchmarkB_ws);
	size += QProblemB_calculateMemorySize(nV);        // qp
	size += DenseMatrix_calculateMemorySize(nV, nV);  // H
	size += 1 * nV * sizeof(real_t);             	  // x
	size += 1 * nV * sizeof(real_t);			      // y

	size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
	size += 1 * 64;                // align once to typical cache line size

	return size;
}

char *OQPbenchmarkB_ws_assignMemory( unsigned int nV, OQPbenchmarkB_ws **mem, void *raw_memory )
{
	// char pointer
	char *c_ptr = (char *)raw_memory;

	// assign structures
	*mem = (OQPbenchmarkB_ws *) c_ptr;
	c_ptr += sizeof(OQPbenchmarkB_ws);

	(*mem)->qp = (QProblemB *) c_ptr;
	c_ptr = QProblemB_assignMemory(nV, &((*mem)->qp), c_ptr);

	(*mem)->H = (DenseMatrix *) c_ptr;
	c_ptr = DenseMatrix_assignMemory(nV, nV, &((*mem)->H), c_ptr);

	// align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
	c_ptr = (char *)s_ptr;

	// assign data
	(*mem)->x = (real_t *) c_ptr;
	c_ptr += nV * sizeof(real_t);

	(*mem)->y = (real_t *) c_ptr;
	c_ptr += nV * sizeof(real_t);

	return c_ptr;
}

OQPbenchmarkB_ws *OQPbenchmarkB_ws_createMemory( unsigned int nV )
{
	OQPbenchmarkB_ws *mem;
    int memory_size = OQPbenchmarkB_ws_calculateMemorySize(nV);
    void *raw_memory_ptr = malloc(memory_size);
    char *ptr_end =  OQPbenchmarkB_ws_assignMemory(nV, &mem, raw_memory_ptr);
    assert((char*)raw_memory_ptr + memory_size >= ptr_end); (void) ptr_end;
    return mem;
}

int OQPinterface_ws_calculateMemorySize( unsigned int nV, unsigned int nC, unsigned int nQP  )
{
	int size = 0;
	size += sizeof(OQPinterface_ws);					  // structure itself
	size += OQPbenchmark_ws_calculateMemorySize(nV, nC);  // qp_ws
	size += OQPbenchmarkB_ws_calculateMemorySize(nV);     // qpB_ws
	size += (nV * nV) * sizeof(real_t);					  // H
	size += (nQP * nV) * sizeof(real_t);				  // g
	size += (nC * nV) * sizeof(real_t);					  // A
	size += (nQP * nV) * sizeof(real_t);				  // lb
	size += (nQP * nV) * sizeof(real_t);				  // ub
	size += (nQP * nC) * sizeof(real_t);				  // lbA
	size += (nQP * nC) * sizeof(real_t);				  // ubA

	size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
	size += 1 * 64;                // align once to typical cache line size

	return size;
}

char *OQPinterface_ws_assignMemory( unsigned int nV, unsigned int nC, unsigned int nQP, OQPinterface_ws **mem, void *raw_memory )
{
	// char pointer
	char *c_ptr = (char *)raw_memory;

	// assign structures
	*mem = (OQPinterface_ws *) c_ptr;
	c_ptr += sizeof(OQPinterface_ws);

	(*mem)->qp_ws = (OQPbenchmark_ws *) c_ptr;
	c_ptr = OQPbenchmark_ws_assignMemory(nV, nC, &((*mem)->qp_ws), c_ptr);

	(*mem)->qpB_ws = (OQPbenchmarkB_ws *) c_ptr;
	c_ptr = OQPbenchmarkB_ws_assignMemory(nV, &((*mem)->qpB_ws), c_ptr);

	// align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
	c_ptr = (char *)s_ptr;

	// assign data
	(*mem)->H = (real_t *) c_ptr;
	c_ptr += (nV * nV) * sizeof(real_t);

	(*mem)->g = (real_t *) c_ptr;
	c_ptr += (nQP * nV) * sizeof(real_t);

	(*mem)->A = (real_t *) c_ptr;
	c_ptr += (nC * nV) * sizeof(real_t);

	(*mem)->lb = (real_t *) c_ptr;
	c_ptr += (nQP * nV) * sizeof(real_t);

	(*mem)->ub = (real_t *) c_ptr;
	c_ptr += (nQP * nV) * sizeof(real_t);

	(*mem)->lbA = (real_t *) c_ptr;
	c_ptr += (nQP * nC) * sizeof(real_t);

	(*mem)->ubA = (real_t *) c_ptr;
	c_ptr += (nQP * nC) * sizeof(real_t);

	return c_ptr;
}

OQPinterface_ws *OQPinterface_ws_createMemory( unsigned int nV, unsigned int nC, unsigned int nQP )
{
	OQPinterface_ws *mem;
    int memory_size = OQPinterface_ws_calculateMemorySize(nV, nC, nQP);
    void *raw_memory_ptr = malloc(memory_size);
    char *ptr_end =  OQPinterface_ws_assignMemory(nV, nC, nQP, &mem, raw_memory_ptr);
    assert((char*)raw_memory_ptr + memory_size >= ptr_end); (void) ptr_end;
    return mem;
}

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
								real_t* maxStationarity, real_t* maxFeasibility, real_t* maxComplementarity,
								OQPbenchmark_ws* work
								)
{
	int k;

	QProblem *qp = work->qp;
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
	real_t *x = work->x;
	real_t *y = work->y;

	/* 4) Prepare matrix objects */
	DenseMatrix *H = work->H;
	DenseMatrix *A = work->A;

	DenseMatrixCON( H, nV, nV, nV, _H );
	DenseMatrixCON( A, nC, nV, nV, _A );

	*maxNWSR = 0;
	*avgNWSR = 0;
	*maxCPUtime = 0.0;
	*avgCPUtime = 0.0;
	*maxStationarity = 0.0;
	*maxFeasibility = 0.0;
	*maxComplementarity = 0.0;

	/*DenseMatrix_print( H );*/

	/* II) SETUP QPROBLEM OBJECT */
	QProblemCON( qp,nV,nC,HST_UNKNOWN );
	QProblem_setOptions( qp,*options );
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
			returnvalue = QProblem_initM( qp, H,gCur,A,lbCur,ubCur,lbACur,ubACur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* hotstart */
			returnvalue = QProblem_hotstart( qp, gCur,lbCur,ubCur,lbACur,ubACur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}

		/* 4) Obtain solution vectors and objective function value */
		QProblem_getPrimalSolution( qp,x );
		QProblem_getDualSolution( qp,y );

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
								real_t* maxStationarity, real_t* maxFeasibility, real_t* maxComplementarity,
								OQPbenchmarkB_ws *work
								)
{
	int k;

	QProblemB *qp = work->qp;
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
	real_t *x = work->x;
	real_t *y = work->y;

	/* 4) Prepare matrix objects */
	DenseMatrix *H = work->H;

	DenseMatrixCON( H, nV, nV, nV, _H );

	*maxNWSR = 0;
	*avgNWSR = 0;
	*maxCPUtime = 0.0;
	*avgCPUtime = 0.0;
	*maxStationarity = 0.0;
	*maxFeasibility = 0.0;
	*maxComplementarity = 0.0;

	/* II) SETUP QPROBLEM OBJECT */
	QProblemBCON( qp,nV,HST_UNKNOWN );
	QProblemB_setOptions( qp,*options );
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
			returnvalue = QProblemB_initM( qp,H,gCur,lbCur,ubCur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}
		else
		{
			/* hotstart */
			returnvalue = QProblemB_hotstart( qp,gCur,lbCur,ubCur, &nWSRcur,&CPUtimeCur );
			if ( ( returnvalue != SUCCESSFUL_RETURN ) && ( returnvalue != RET_MAX_NWSR_REACHED ) )
				return THROWERROR( returnvalue );
		}

		/* 4) Obtain solution vectors and objective function value ... */
		QProblemB_getPrimalSolution( qp,x );
		QProblemB_getDualSolution( qp,y );

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
								real_t* maxStationarity, real_t* maxFeasibility, real_t* maxComplementarity,
								OQPinterface_ws* work
								)
{
	int nQP=0, nV=0, nC=0, nEC=0;

	real_t *H = work->H;
	real_t *g = work->g;
	real_t *A = work->A;
	real_t *lb = work->lb;
	real_t *ub = work->ub;
	real_t *lbA = work->lbA;
	real_t *ubA = work->ubA;

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
											maxStationarity,maxFeasibility,maxComplementarity,
											work->qp_ws
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
											maxStationarity,maxFeasibility,maxComplementarity,
											work->qpB_ws
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

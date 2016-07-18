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
 *	\file interfaces/matlab/qpOASES_e_matlab_utils.c
 *	\author Hans Joachim Ferreau, Alexander Buchner
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Collects utility functions for Interface to Matlab(R) that
 *	enables to call qpOASES as a MEX function.
 *
 */


void QPInstanceCON(	QPInstance* _THIS,
					int _nV,
					int _nC,
					HessianType _hessianType,
					BooleanType _isSimplyBounded
					)
{
	if ( _nV > NVMAX )
	{
		myMexErrMsgTxt( "ERROR (qpOASES_e): Too many primal variables (try increasing QPOASES_NVMAX)!" );
		return;
	}

	if ( _nC > NCMAX )
	{
		myMexErrMsgTxt( "ERROR (qpOASES_e): Too many constraints (try increasing QPOASES_NCMAX)!" );
		return;
	}

	_THIS->handle = QPInstance_nexthandle;
	
	if ( _nC > 0 )
		_THIS->isSimplyBounded = BT_FALSE;
	else
		_THIS->isSimplyBounded = _isSimplyBounded;
	
	if ( _THIS->isSimplyBounded == BT_FALSE )
		QProblemCON( &(_THIS->sqp),_nV,_nC,_hessianType );
	else
		QProblemBCON( &(_THIS->qpb),_nV,_hessianType );
}	


int QPInstance_getNV( QPInstance* _THIS )
{
    if ( _THIS->isSimplyBounded == BT_FALSE )
        return QProblem_getNV( &(_THIS->sqp) );
	else
		return QProblemB_getNV( &(_THIS->qpb) );
}


int QPInstance_getNC( QPInstance* _THIS )
{
    if ( _THIS->isSimplyBounded == BT_FALSE )
        return QProblem_getNC( &(_THIS->sqp) );
	else
		return 0;
}



/*
 *	m x I s S c a l a r
 */
bool mxIsScalar( const mxArray *pm )
{
	if ( ( mxGetM(pm) == 1 ) && ( mxGetN(pm) == 1 ) )
		return true;
	else
		return false;
}



/*
 *	a l l o c a t e Q P r o b l e m I n s t a n c e
 */
int allocateQPInstance(	int nV, int nC, HessianType hessianType,
						BooleanType isSimplyBounded, const Options* options
						)
{
	QPInstance* inst = 0;

	QPInstance_nexthandle = (QPInstance_nexthandle+1) % MAX_NUM_QPINSTANCES;

	inst = &(QPInstances[QPInstance_nexthandle]);
	QPInstanceCON( inst, nV,nC,hessianType, isSimplyBounded );

	if ( options != 0 )
	{
		if ( isSimplyBounded == BT_FALSE )
			QProblem_setOptions( &(inst->sqp),*options );
		else
			QProblemB_setOptions( &(inst->qpb),*options );
	}
	
	return inst->handle;
}


/*
 *  g e t Q P r o b l e m I n s t a n c e
 */
QPInstance* getQPInstance( int handle )
{
	if ( handle < MAX_NUM_QPINSTANCES )
		return &(QPInstances[handle]);
	else
		return 0;
}



/*
 *	s m a r t D i m e n s i o n C h e c k
 */
returnValue smartDimensionCheck(	real_t** input, unsigned int m, unsigned int n, BooleanType emptyAllowed,
									const mxArray* prhs[], int idx
									)
{
	char msg[QPOASES_MAX_STRING_LENGTH];

	/* If index is negative, the input does not exist. */
	if ( idx < 0 )
	{
		*input = 0;
		return SUCCESSFUL_RETURN;
	}

	/* Otherwise the input has been passed by the user. */
	if ( mxIsEmpty( prhs[ idx ] ) )
	{
		/* input is empty */
		if ( ( emptyAllowed == BT_TRUE ) || ( idx == 0 ) ) /* idx==0 used for auxInput */
		{
			*input = 0;
			return SUCCESSFUL_RETURN;
		}
		else
		{
			if ( idx > 0 )
				snprintf(msg, QPOASES_MAX_STRING_LENGTH, "ERROR (qpOASES_e): Empty argument %d not allowed!", idx+1);
			myMexErrMsgTxt( msg );
			return RET_INVALID_ARGUMENTS;
		}
	}
	else
	{
		/* input is non-empty */
        if ( mxIsSparse( prhs[ idx ] ) == 0 )
        {
            if ( ( mxGetM( prhs[ idx ] ) == m ) && ( mxGetN( prhs[ idx ] ) == n ) )
            {
                *input = (real_t*) mxGetPr( prhs[ idx ] );
                return SUCCESSFUL_RETURN;
            }
            else
            {
				if ( idx > 0 )
					snprintf(msg, QPOASES_MAX_STRING_LENGTH, "ERROR (qpOASES_e): Input dimension mismatch for argument %d ([%ld,%ld] ~= [%d,%d]).",
							 idx+1, (long int)mxGetM(prhs[idx]), (long int)mxGetN(prhs[idx]), m, n);
				else /* idx==0 used for auxInput */
					snprintf(msg, QPOASES_MAX_STRING_LENGTH, "ERROR (qpOASES_e): Input dimension mismatch for some auxInput entry ([%ld,%ld] ~= [%d,%d]).",
							 (long int)mxGetM(prhs[idx]), (long int)mxGetN(prhs[idx]), m, n);
                myMexErrMsgTxt( msg );
                return RET_INVALID_ARGUMENTS;
            }
        }
        else
        {
			if ( idx > 0 )
				snprintf(msg, QPOASES_MAX_STRING_LENGTH, "ERROR (qpOASES_e): Vector argument %d must not be in sparse format!", idx+1);
			else /* idx==0 used for auxInput */
				snprintf(msg, QPOASES_MAX_STRING_LENGTH, "ERROR (qpOASES_e): auxInput entries must not be in sparse format!" );
			myMexErrMsgTxt( msg );
			return RET_INVALID_ARGUMENTS;
        }
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	c o n t a i n s N a N
 */
BooleanType containsNaN( const real_t* const data, unsigned int dim )
{
	unsigned int i;

	if ( data == 0 )
		return BT_FALSE;

	for ( i = 0; i < dim; ++i )
		if ( mxIsNaN(data[i]) == 1 )
			return BT_TRUE;

	return BT_FALSE;
}


/*
 *	c o n t a i n s I n f
 */
BooleanType containsInf( const real_t* const data, unsigned int dim )
{
	unsigned int i;

	if ( data == 0 )
		return BT_FALSE;

	for ( i = 0; i < dim; ++i )
		if ( mxIsInf(data[i]) == 1 )
			return BT_TRUE;

	return BT_FALSE;
}


/*
 *	c o n t a i n s N a N o r I n f
 */
BooleanType containsNaNorInf(	const mxArray* prhs[], int rhs_index,
								BooleanType mayContainInf
								)
{
	unsigned int dim;
	char msg[QPOASES_MAX_STRING_LENGTH];

	if ( rhs_index < 0 )
		return BT_FALSE;

	/* overwrite dim for sparse matrices */
	if (mxIsSparse(prhs[rhs_index]) == 1)
		dim = (unsigned int)mxGetNzmax(prhs[rhs_index]);
	else
		dim = mxGetM(prhs[rhs_index]) * mxGetN(prhs[rhs_index]);

	if (containsNaN((real_t*) mxGetPr(prhs[rhs_index]), dim) == BT_TRUE) {
		snprintf(msg, QPOASES_MAX_STRING_LENGTH,
				"ERROR (qpOASES_e): Argument %d contains 'NaN' !", rhs_index + 1);
		myMexErrMsgTxt(msg);
		return BT_TRUE;
	}

	if ( mayContainInf == BT_FALSE ) {
		if (containsInf((real_t*) mxGetPr(prhs[rhs_index]), dim) == BT_TRUE) {
			snprintf(msg, QPOASES_MAX_STRING_LENGTH,
					"ERROR (qpOASES_e): Argument %d contains 'Inf' !",
					rhs_index + 1);
			myMexErrMsgTxt(msg);
			return BT_TRUE;
		}
	}

	return BT_FALSE;
}


/*
 *	c o n v e r t F o r t r a n T o C
 */
returnValue convertFortranToC( const real_t* const M_for, int nV, int nC, real_t* const M )
{
	int i,j;

	if ( ( M_for == 0 ) || ( M == 0 ) )
		return RET_INVALID_ARGUMENTS;

	if ( ( nV < 0 ) || ( nC < 0 ) )
		return RET_INVALID_ARGUMENTS;

	for ( i=0; i<nC; ++i )
		for ( j=0; j<nV; ++j )
			M[i*nV + j] = M_for[j*nC + i];

	return SUCCESSFUL_RETURN;
}


/*
 *	h a s O p t i o n s V a l u e
 */
BooleanType hasOptionsValue( const mxArray* optionsPtr, const char* const optionString, double** optionValue )
{
	mxArray* optionName = mxGetField( optionsPtr,0,optionString );

	if ( optionName == 0 )
	{
		char msg[QPOASES_MAX_STRING_LENGTH];
		snprintf(msg, QPOASES_MAX_STRING_LENGTH, "Option struct does not contain entry '%s', using default value instead!", optionString );
		mexWarnMsgTxt( msg );
		return BT_FALSE;
	}

	if ( ( mxIsEmpty(optionName) == false ) && ( mxIsScalar( optionName ) == true ) )
	{
		*optionValue = mxGetPr( optionName );
		return BT_TRUE;
	}
	else
	{
		char msg[QPOASES_MAX_STRING_LENGTH];
		snprintf(msg, QPOASES_MAX_STRING_LENGTH, "Option '%s' is not a scalar, using default value instead!", optionString );
		mexWarnMsgTxt( msg );
		return BT_FALSE;
	}
}


/*
 *	s e t u p O p t i o n s
 */
returnValue setupOptions( Options* options, const mxArray* optionsPtr, int* nWSRin, real_t* maxCpuTime )
{
	double* optionValue;
	int optionValueInt;

	/* Check for correct number of option entries;
	 * may occur, e.g., if user types options.<misspelledName> = <someValue>; */
	if ( mxGetNumberOfFields(optionsPtr) != 31 )
		mexWarnMsgTxt( "Options might be set incorrectly as struct has wrong number of entries!\n         Type 'help qpOASES_options' for further information." );


	if ( hasOptionsValue( optionsPtr,"maxIter",&optionValue ) == BT_TRUE )
		if ( *optionValue >= 0.0 )
			*nWSRin = (int)*optionValue;

	if ( hasOptionsValue( optionsPtr,"maxCpuTime",&optionValue ) == BT_TRUE )
		if ( *optionValue >= 0.0 )
			*maxCpuTime = *optionValue;

	if ( hasOptionsValue( optionsPtr,"printLevel",&optionValue ) == BT_TRUE )
	{
        #ifdef __SUPPRESSANYOUTPUT__
        options->printLevel = PL_NONE;
        #else
		optionValueInt = (int)*optionValue;
		options->printLevel = (REFER_NAMESPACE_QPOASES PrintLevel)optionValueInt;
        if ( options->printLevel < PL_DEBUG_ITER )
            options->printLevel = PL_DEBUG_ITER;
        if ( options->printLevel > PL_HIGH )
            options->printLevel = PL_HIGH;       
        #endif
	}

	if ( hasOptionsValue( optionsPtr,"enableRamping",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		options->enableRamping = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
	}

	if ( hasOptionsValue( optionsPtr,"enableFarBounds",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		options->enableFarBounds = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
	}

	if ( hasOptionsValue( optionsPtr,"enableFlippingBounds",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		options->enableFlippingBounds = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
	}

	if ( hasOptionsValue( optionsPtr,"enableRegularisation",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		options->enableRegularisation = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
	}

	if ( hasOptionsValue( optionsPtr,"enableFullLITests",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		options->enableFullLITests = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
	}

	if ( hasOptionsValue( optionsPtr,"enableNZCTests",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		options->enableNZCTests = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
	}

	if ( hasOptionsValue( optionsPtr,"enableDriftCorrection",&optionValue ) == BT_TRUE )
		options->enableDriftCorrection = (int)*optionValue;

	if ( hasOptionsValue( optionsPtr,"enableCholeskyRefactorisation",&optionValue ) == BT_TRUE )
		options->enableCholeskyRefactorisation = (int)*optionValue;

	if ( hasOptionsValue( optionsPtr,"enableEqualities",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		options->enableEqualities = (REFER_NAMESPACE_QPOASES BooleanType)optionValueInt;
	}


	if ( hasOptionsValue( optionsPtr,"terminationTolerance",&optionValue ) == BT_TRUE )
		options->terminationTolerance = *optionValue;

	if ( hasOptionsValue( optionsPtr,"boundTolerance",&optionValue ) == BT_TRUE )
		options->boundTolerance = *optionValue;

	if ( hasOptionsValue( optionsPtr,"boundRelaxation",&optionValue ) == BT_TRUE )
		options->boundRelaxation = *optionValue;

	if ( hasOptionsValue( optionsPtr,"epsNum",&optionValue ) == BT_TRUE )
		options->epsNum = *optionValue;

	if ( hasOptionsValue( optionsPtr,"epsDen",&optionValue ) == BT_TRUE )
		options->epsDen = *optionValue;

	if ( hasOptionsValue( optionsPtr,"maxPrimalJump",&optionValue ) == BT_TRUE )
		options->maxPrimalJump = *optionValue;

	if ( hasOptionsValue( optionsPtr,"maxDualJump",&optionValue ) == BT_TRUE )
		options->maxDualJump = *optionValue;


	if ( hasOptionsValue( optionsPtr,"initialRamping",&optionValue ) == BT_TRUE )
		options->initialRamping = *optionValue;

	if ( hasOptionsValue( optionsPtr,"finalRamping",&optionValue ) == BT_TRUE )
		options->finalRamping = *optionValue;

	if ( hasOptionsValue( optionsPtr,"initialFarBounds",&optionValue ) == BT_TRUE )
		options->initialFarBounds = *optionValue;

	if ( hasOptionsValue( optionsPtr,"growFarBounds",&optionValue ) == BT_TRUE )
		options->growFarBounds = *optionValue;

	if ( hasOptionsValue( optionsPtr,"initialStatusBounds",&optionValue ) == BT_TRUE )
	{
		optionValueInt = (int)*optionValue;
		if ( optionValueInt < -1 ) 
			optionValueInt = -1;
		if ( optionValueInt > 1 ) 
			optionValueInt = 1;
		options->initialStatusBounds = (REFER_NAMESPACE_QPOASES SubjectToStatus)optionValueInt;
	}

	if ( hasOptionsValue( optionsPtr,"epsFlipping",&optionValue ) == BT_TRUE )
		options->epsFlipping = *optionValue;

	if ( hasOptionsValue( optionsPtr,"numRegularisationSteps",&optionValue ) == BT_TRUE )
		options->numRegularisationSteps = (int)*optionValue;

	if ( hasOptionsValue( optionsPtr,"epsRegularisation",&optionValue ) == BT_TRUE )
		options->epsRegularisation = *optionValue;

	if ( hasOptionsValue( optionsPtr,"numRefinementSteps",&optionValue ) == BT_TRUE )
		options->numRefinementSteps = (int)*optionValue;

	if ( hasOptionsValue( optionsPtr,"epsIterRef",&optionValue ) == BT_TRUE )
		options->epsIterRef = *optionValue;

	if ( hasOptionsValue( optionsPtr,"epsLITests",&optionValue ) == BT_TRUE )
		options->epsLITests = *optionValue;

	if ( hasOptionsValue( optionsPtr,"epsNZCTests",&optionValue ) == BT_TRUE )
		options->epsNZCTests = *optionValue;

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p A u x i l i a r y I n p u t s
 */
returnValue setupAuxiliaryInputs(	const mxArray* auxInput, unsigned int nV, unsigned int nC,
									HessianType* hessianType, double** x0, double** guessedBounds, double** guessedConstraints, double** R
									)
{
	mxArray* curField = 0;
	double* hessianTypeTmp;
	int hessianTypeInt;

	/* hessianType */
	curField = mxGetField( auxInput,0,"hessianType" );
	if ( curField == NULL )
		mexWarnMsgTxt( "auxInput struct does not contain entry 'hessianType'!\n         Type 'help qpOASES_auxInput' for further information." );
	else
	{
		if ( mxIsEmpty(curField) == true )
		{
			*hessianType = HST_UNKNOWN;
		}
		else
		{
			if ( mxIsScalar(curField) == false )
				return RET_INVALID_ARGUMENTS;

			hessianTypeTmp = mxGetPr(curField);
			hessianTypeInt = (int)*hessianTypeTmp;
			if ( hessianTypeInt < 0 ) 
				hessianTypeInt = 6; /* == HST_UNKNOWN */
			if ( hessianTypeInt > 5 ) 
				hessianTypeInt = 6; /* == HST_UNKNOWN */
			*hessianType = (REFER_NAMESPACE_QPOASES HessianType)hessianTypeInt;
		}
	}

	/* x0 */
	curField = mxGetField( auxInput,0,"x0" );
	if ( curField == NULL )
		mexWarnMsgTxt( "auxInput struct does not contain entry 'x0'!\n         Type 'help qpOASES_auxInput' for further information." );
	else
	{
		*x0 = mxGetPr(curField);
		if ( smartDimensionCheck( x0,nV,1, BT_TRUE,((const mxArray**)&curField),0 ) != SUCCESSFUL_RETURN )
			return RET_INVALID_ARGUMENTS;
	}

	/* guessedWorkingSetB */
	curField = mxGetField( auxInput,0,"guessedWorkingSetB" );
	if ( curField == NULL )
		mexWarnMsgTxt( "auxInput struct does not contain entry 'guessedWorkingSetB'!\n         Type 'help qpOASES_auxInput' for further information." );
	else
	{
		*guessedBounds = mxGetPr(curField);
		if ( smartDimensionCheck( guessedBounds,nV,1, BT_TRUE,((const mxArray**)&curField),0 ) != SUCCESSFUL_RETURN )
			return RET_INVALID_ARGUMENTS;
	}

	/* guessedWorkingSetC */
	curField = mxGetField( auxInput,0,"guessedWorkingSetC" );
	if ( curField == NULL )
		mexWarnMsgTxt( "auxInput struct does not contain entry 'guessedWorkingSetC'!\n         Type 'help qpOASES_auxInput' for further information." );
	else
	{
		*guessedConstraints = mxGetPr(curField);
		if ( smartDimensionCheck( guessedConstraints,nC,1, BT_TRUE,((const mxArray**)&curField),0 ) != SUCCESSFUL_RETURN )
			return RET_INVALID_ARGUMENTS;
	}

	/* R */
	curField = mxGetField( auxInput,0,"R" );
	if ( curField == NULL )
		mexWarnMsgTxt( "auxInput struct does not contain entry 'R'!\n         Type 'help qpOASES_auxInput' for further information." );
	else
	{
		*R = mxGetPr(curField);
		if ( smartDimensionCheck( R,nV,nV, BT_TRUE,((const mxArray**)&curField),0 ) != SUCCESSFUL_RETURN )
			return RET_INVALID_ARGUMENTS;
	}

	return SUCCESSFUL_RETURN;
}



/*
 *	a l l o c a t e O u t p u t s
 */
returnValue allocateOutputs(	int nlhs, mxArray* plhs[], int nV, int nC, int nP, int handle
								)
{
	/* Create output vectors and assign pointers to them. */
	int curIdx = 0;

	/* handle */
	if ( handle >= 0 )
		plhs[curIdx++] = mxCreateDoubleMatrix( 1, 1, mxREAL );

	/* x */
	plhs[curIdx++] = mxCreateDoubleMatrix( nV, nP, mxREAL );

	if ( nlhs > curIdx )
	{
		/* fval */
		plhs[curIdx++] = mxCreateDoubleMatrix( 1, nP, mxREAL );

		if ( nlhs > curIdx )
		{
			/* exitflag */
			plhs[curIdx++] = mxCreateDoubleMatrix( 1, nP, mxREAL );

			if ( nlhs > curIdx )
			{
				/* iter */
				plhs[curIdx++] = mxCreateDoubleMatrix( 1, nP, mxREAL );

				if ( nlhs > curIdx )
				{
					/* lambda */
					plhs[curIdx++] = mxCreateDoubleMatrix( nV+nC, nP, mxREAL );

					if ( nlhs > curIdx )
					{
						/* setup auxiliary output struct */
						mxArray* auxOutput = mxCreateStructMatrix( 1,1,0,0 );
						int curFieldNum;
						
						/* working set */
						curFieldNum = mxAddField( auxOutput,"workingSetB" );
						if ( curFieldNum >= 0 )
							mxSetFieldByNumber( auxOutput,0,curFieldNum,mxCreateDoubleMatrix( nV, nP, mxREAL ) );

						curFieldNum = mxAddField( auxOutput,"workingSetC" );
						if ( curFieldNum >= 0 )
							mxSetFieldByNumber( auxOutput,0,curFieldNum,mxCreateDoubleMatrix( nC, nP, mxREAL ) );

						curFieldNum = mxAddField( auxOutput,"cpuTime" );
						if ( curFieldNum >= 0 )
							mxSetFieldByNumber( auxOutput,0,curFieldNum,mxCreateDoubleMatrix( 1, nP, mxREAL ) );

						plhs[curIdx] = auxOutput;
					}
				}
			}
		}
	}
	
	return SUCCESSFUL_RETURN;
}


/*
 *	o b t a i n O u t p u t s
 */
returnValue obtainOutputs(	int k, QProblem* qp, returnValue returnvalue, int _nWSRout, double _cpuTime,
							int nlhs, mxArray* plhs[], int nV, int nC, int handle
							)
{
	double *x=0, *obj=0, *status=0, *nWSRout=0, *y=0, *workingSetB=0, *workingSetC=0, *cpuTime;
	mxArray *auxOutput=0, *curField=0;

	/* Create output vectors and assign pointers to them. */
	int curIdx = 0;

	/* handle */
	if ( handle >= 0 )
		plhs[curIdx++] = mxCreateDoubleScalar( handle );

	/* x */
	x = mxGetPr( plhs[curIdx++] );
	QProblem_getPrimalSolution( qp,&(x[k*nV]) );

	if ( nlhs > curIdx )
	{
		/* fval */
		obj = mxGetPr( plhs[curIdx++] );
		obj[k] = QProblem_getObjVal( qp );

		if ( nlhs > curIdx )
		{
			/* exitflag */
			status = mxGetPr( plhs[curIdx++] );
			status[k] = (double)qpOASES_getSimpleStatus( returnvalue,0 );

			if ( nlhs > curIdx )
			{
				/* iter */
				nWSRout = mxGetPr( plhs[curIdx++] );
				nWSRout[k] = (double) _nWSRout;

				if ( nlhs > curIdx )
				{
					/* lambda */
					y = mxGetPr( plhs[curIdx++] );
					QProblem_getDualSolution( qp,&(y[k*(nV+nC)]) );

					/* auxOutput */
					if ( nlhs > curIdx )
					{
						auxOutput = plhs[curIdx];
						curField = 0;

						/* working set bounds */
						if ( nV > 0 )
						{
							curField = mxGetField( auxOutput,0,"workingSetB" );
							workingSetB = mxGetPr(curField);
							QProblem_getWorkingSetBounds( qp,&(workingSetB[k*nV]) );
						}

						/* working set constraints */
						if ( nC > 0 )
						{
							curField = mxGetField( auxOutput,0,"workingSetC" );
							workingSetC = mxGetPr(curField);
							QProblem_getWorkingSetConstraints( qp,&(workingSetC[k*nC]) );
						}

						/* cpu time */
						curField = mxGetField( auxOutput,0,"cpuTime" );
						cpuTime = mxGetPr(curField);
						cpuTime[0] = (double) _cpuTime;
					}
				}
			}
		}
	}
	
	return SUCCESSFUL_RETURN;
}


/*
 *	o b t a i n O u t p u t s S B
 */
returnValue obtainOutputsSB(	int k, QProblemB* qp, returnValue returnvalue, int _nWSRout, double _cpuTime,
								int nlhs, mxArray* plhs[], int nV, int handle
								)
{
	double *x=0, *obj=0, *status=0, *nWSRout=0, *y=0, *workingSetB=0, *workingSetC=0, *cpuTime;
	mxArray *auxOutput=0, *curField=0;

	/* Create output vectors and assign pointers to them. */
	int curIdx = 0;

	/* handle */
	if ( handle >= 0 )
		plhs[curIdx++] = mxCreateDoubleScalar( handle );

	/* x */
	x = mxGetPr( plhs[curIdx++] );
	QProblemB_getPrimalSolution( qp,&(x[k*nV]) );

	if ( nlhs > curIdx )
	{
		/* fval */
		obj = mxGetPr( plhs[curIdx++] );
		obj[k] = QProblemB_getObjVal( qp );

		if ( nlhs > curIdx )
		{
			/* exitflag */
			status = mxGetPr( plhs[curIdx++] );
			status[k] = (double)qpOASES_getSimpleStatus( returnvalue,0 );

			if ( nlhs > curIdx )
			{
				/* iter */
				nWSRout = mxGetPr( plhs[curIdx++] );
				nWSRout[k] = (double) _nWSRout;

				if ( nlhs > curIdx )
				{
					/* lambda */
					y = mxGetPr( plhs[curIdx++] );
					QProblemB_getDualSolution( qp,&(y[k*nV]) );

					/* auxOutput */
					if ( nlhs > curIdx )
					{
						auxOutput = plhs[curIdx];
						curField = 0;

						/* working set bounds */
						if ( nV > 0 )
						{
							curField = mxGetField( auxOutput,0,"workingSetB" );
							workingSetB = mxGetPr(curField);
							QProblemB_getWorkingSetBounds( qp,&(workingSetB[k*nV]) );
						}

						/* cpu time */
						curField = mxGetField( auxOutput,0,"cpuTime" );
						cpuTime = mxGetPr(curField);
						cpuTime[0] = (double) _cpuTime;
					}
				}
			}
		}
	}
	
	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p H e s s i a n M a t r i x
 */
returnValue setupHessianMatrix(	const mxArray* prhsH, int nV,
								DenseMatrix* H
								)
{
	real_t *H_for = 0;
	/*static real_t H_mem[NVMAX*NVMAX];*/

	if ( prhsH == 0 )
		return SUCCESSFUL_RETURN;

	if ( mxIsSparse( prhsH ) != 0 )
	{
		myMexErrMsgTxt( "ERROR (qpOASES_e): Cannot handle sparse matrices!" );
		return RET_INVALID_ARGUMENTS;
	}
	else
	{
		/* make a deep-copy in order to avoid modifying input data when regularising */
		H_for = (real_t*) mxGetPr( prhsH );
		/*qpOASES_printM( H_for,nV,nV );*/
		/*convertFortranToC( H_for, nV,nV,H_mem );  not neccessary as H is symmetric! */
		DenseMatrixCON( H, nV,nV,nV, H_for );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p C o n s t r a i n t M a t r i x
 */
returnValue setupConstraintMatrix(	const mxArray* prhsA, int nV, int nC,
									DenseMatrix* A
									)
{
	real_t *A_for = 0;
	static real_t A_mem[NCMAX*NVMAX];

	if ( prhsA == 0 )
		return SUCCESSFUL_RETURN;

	if ( mxIsSparse( prhsA ) != 0 )
	{
		myMexErrMsgTxt( "ERROR (qpOASES_e): Cannot handle sparse matrices!" );
		return RET_INVALID_ARGUMENTS;
	}
	else
	{
		/* make a deep-copy in order to avoid modifying input data when regularising */
		A_for = (real_t*) mxGetPr( prhsA );
		convertFortranToC( A_for, nV,nC,A_mem );
		DenseMatrixCON( A, nC,nV,nV, A_mem );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	end of file
 */

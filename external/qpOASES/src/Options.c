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
 *	\file src/Options.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Implementation of the Options class designed to manage working sets of
 *	constraints and bounds within a QProblem.
 */


#include <qpOASES_e/Options.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	O p t i o n s
 */
void OptionsCON(	Options* _THIS
					)
{
	#ifdef __CODE_GENERATION__
	Options_setToFast( _THIS );
	#else
	Options_setToDefault( _THIS );
	#endif /* __CODE_GENERATION__ */
}


/*
 *	c o p y
 */
void OptionsCPY(	Options* FROM,
					Options* TO
					)
{
	TO->printLevel                    =  FROM->printLevel;

	TO->enableRamping                 =  FROM->enableRamping;
	TO->enableFarBounds               =  FROM->enableFarBounds;
	TO->enableFlippingBounds          =  FROM->enableFlippingBounds;
	TO->enableRegularisation          =  FROM->enableRegularisation;
	TO->enableFullLITests             =  FROM->enableFullLITests;
	TO->enableNZCTests                =  FROM->enableNZCTests;
	TO->enableDriftCorrection         =  FROM->enableDriftCorrection;
	TO->enableCholeskyRefactorisation =  FROM->enableCholeskyRefactorisation;
	TO->enableEqualities              =  FROM->enableEqualities;

	TO->terminationTolerance          =  FROM->terminationTolerance;
	TO->boundTolerance                =  FROM->boundTolerance;
	TO->boundRelaxation               =  FROM->boundRelaxation;
	TO->epsNum                        =  FROM->epsNum;
	TO->epsDen                        =  FROM->epsDen;
	TO->maxPrimalJump                 =  FROM->maxPrimalJump;
	TO->maxDualJump                   =  FROM->maxDualJump;

	TO->initialRamping                =  FROM->initialRamping;
	TO->finalRamping                  =  FROM->finalRamping;
	TO->initialFarBounds              =  FROM->initialFarBounds;
	TO->growFarBounds                 =  FROM->growFarBounds;
 	TO->initialStatusBounds           =  FROM->initialStatusBounds;
	TO->epsFlipping                   =  FROM->epsFlipping;
	TO->numRegularisationSteps        =  FROM->numRegularisationSteps;
	TO->epsRegularisation             =  FROM->epsRegularisation;
	TO->numRefinementSteps            =  FROM->numRefinementSteps;
	TO->epsIterRef                    =  FROM->epsIterRef;
	TO->epsLITests                    =  FROM->epsLITests;
	TO->epsNZCTests                   =  FROM->epsNZCTests;

	TO->enableDropInfeasibles         =  FROM->enableDropInfeasibles;
    TO->dropBoundPriority             =  FROM->dropBoundPriority;
    TO->dropEqConPriority             =  FROM->dropEqConPriority;
    TO->dropIneqConPriority           =  FROM->dropIneqConPriority;
}



/*
 *	s e t T o D e f a u l t
 */
returnValue Options_setToDefault(	Options* _THIS
									)
{
	_THIS->printLevel = PL_MEDIUM;
	#ifdef __DEBUG__
	_THIS->printLevel = PL_HIGH;
	#endif
	#ifdef __SUPPRESSANYOUTPUT__
	_THIS->printLevel = PL_NONE;
	#endif
	#ifdef __CODE_GENERATION__
	_THIS->printLevel = QPOASES_PRINTLEVEL;
	#endif

	_THIS->enableRamping                 =  BT_TRUE;
	_THIS->enableFarBounds               =  BT_TRUE;
	_THIS->enableFlippingBounds          =  BT_TRUE;
	_THIS->enableRegularisation          =  BT_FALSE;
	_THIS->enableFullLITests             =  BT_FALSE;
	_THIS->enableNZCTests                =  BT_TRUE;
	_THIS->enableDriftCorrection         =  1;
	_THIS->enableCholeskyRefactorisation =  0;
	_THIS->enableEqualities              =  BT_FALSE;

	#ifdef __USE_SINGLE_PRECISION__
	_THIS->terminationTolerance          =  1.0e2 * QPOASES_EPS;
	_THIS->boundTolerance                =  1.0e2 * QPOASES_EPS;
	#else
	_THIS->terminationTolerance          =  5.0e6 * QPOASES_EPS;
	_THIS->boundTolerance                =  1.0e6 * QPOASES_EPS;
	#endif
	_THIS->boundRelaxation               =  1.0e4;
	#ifdef __USE_SINGLE_PRECISION__
	_THIS->epsNum                        = -1.0e2 * QPOASES_EPS;
	_THIS->epsDen                        =  1.0e2 * QPOASES_EPS;
	#else
	_THIS->epsNum                        = -1.0e3 * QPOASES_EPS;
	_THIS->epsDen                        =  1.0e3 * QPOASES_EPS;
	#endif
	_THIS->maxPrimalJump                 =  1.0e8;
	_THIS->maxDualJump                   =  1.0e8;

	_THIS->initialRamping                =  0.5;
	_THIS->finalRamping                  =  1.0;
	_THIS->initialFarBounds              =  1.0e6;
	_THIS->growFarBounds                 =  1.0e3;
 	_THIS->initialStatusBounds           =  ST_LOWER;
	#ifdef __USE_SINGLE_PRECISION__
	_THIS->epsFlipping                   =  5.0e1 * QPOASES_EPS;
	#else
	_THIS->epsFlipping                   =  1.0e3 * QPOASES_EPS;
	#endif
	_THIS->numRegularisationSteps        =  0;
	#ifdef __USE_SINGLE_PRECISION__
	_THIS->epsRegularisation             =  2.0e1 * QPOASES_EPS;
	_THIS->numRefinementSteps            =  2;
	#else
	_THIS->epsRegularisation             =  1.0e3 * QPOASES_EPS;
	_THIS->numRefinementSteps            =  1;
	#endif
	_THIS->epsIterRef                    =  1.0e2 * QPOASES_EPS;
	#ifdef __USE_SINGLE_PRECISION__
	_THIS->epsLITests                    =  5.0e1 * QPOASES_EPS;
	_THIS->epsNZCTests                   =  1.0e2 * QPOASES_EPS;
	#else
	_THIS->epsLITests                    =  1.0e5 * QPOASES_EPS;
	_THIS->epsNZCTests                   =  3.0e3 * QPOASES_EPS;
	#endif

	_THIS->enableDropInfeasibles         =  BT_FALSE;
    _THIS->dropBoundPriority             =  1;
    _THIS->dropEqConPriority             =  1;
    _THIS->dropIneqConPriority           =  1;
    
	return SUCCESSFUL_RETURN;
}


/*
 *	s e t T o R e l i a b l e
 */
returnValue Options_setToReliable(	Options* _THIS
									)
{
	Options_setToDefault( _THIS );

	_THIS->enableFullLITests             =  BT_TRUE;
	_THIS->enableCholeskyRefactorisation =  1;

	#ifdef __USE_SINGLE_PRECISION__
	_THIS->numRefinementSteps            =  3;
	#else
 	_THIS->numRefinementSteps            =  2;
	#endif /*__USE_SINGLE_PRECISION__ */

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t T o M P C
 */
returnValue Options_setToMPC(	Options* _THIS
								)
{
	Options_setToDefault( _THIS );

	_THIS->enableRamping                 =  BT_FALSE;
	_THIS->enableFarBounds               =  BT_TRUE;
	_THIS->enableFlippingBounds          =  BT_FALSE;
	_THIS->enableRegularisation          =  BT_TRUE;
	_THIS->enableNZCTests                =  BT_FALSE;
	_THIS->enableDriftCorrection         =  0;
	#ifdef __USE_SINGLE_PRECISION__
	_THIS->enableEqualities              =  BT_FALSE;
	#else
	_THIS->enableEqualities              =  BT_TRUE;
	#endif

	#ifdef __USE_SINGLE_PRECISION__
	_THIS->terminationTolerance          =  1.0e3 * QPOASES_EPS;
	#else
	_THIS->terminationTolerance          =  1.0e9 * QPOASES_EPS;
	#endif

	_THIS->initialStatusBounds           =  ST_INACTIVE;
	_THIS->numRegularisationSteps        =  1;
	#ifdef __USE_SINGLE_PRECISION__
	_THIS->numRefinementSteps            =  1;
	#else
	_THIS->numRefinementSteps            =  0;
	#endif

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t T o F a s t
 */
returnValue Options_setToFast(	Options* _THIS
								)
{
	return Options_setToMPC( _THIS );
}



/*
 *	e n s u r e C o n s i s t e n c y
 */
returnValue Options_ensureConsistency(	Options* _THIS
										)
{
	BooleanType needToAdjust = BT_FALSE;

	/* flipping bounds require far bounds */
    /* (ckirches) Removed this as per filter's trust region
	if( enableFlippingBounds == BT_TRUE )
		enableFarBounds = BT_TRUE;
    */
		
	if( _THIS->enableDriftCorrection < 0 )
	{
		_THIS->enableDriftCorrection = 0;
		needToAdjust = BT_TRUE;
	}
	
	if( _THIS->enableCholeskyRefactorisation < 0 )
	{
		_THIS->enableCholeskyRefactorisation = 0;
		needToAdjust = BT_TRUE;
	}


	if ( _THIS->terminationTolerance <= 0.0 )
	{
		_THIS->terminationTolerance = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->epsIterRef <= 0.0 )
	{
		_THIS->epsIterRef = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->epsRegularisation <= 0.0 )
	{
		_THIS->epsRegularisation = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->boundTolerance <= 0.0 )
	{
		_THIS->boundTolerance = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->boundRelaxation <= 0.0 )
	{
		_THIS->boundRelaxation = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}
	
	if ( _THIS->maxPrimalJump <= 0.0 )
	{
		_THIS->maxPrimalJump = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->maxDualJump <= 0.0 )
	{
		_THIS->maxDualJump = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}


	if ( _THIS->initialRamping < 0.0 )
	{
		_THIS->initialRamping = 0.0;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->finalRamping < 0.0 )
	{
		_THIS->finalRamping = 0.0;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->initialFarBounds <= _THIS->boundRelaxation )
	{
		_THIS->initialFarBounds = _THIS->boundRelaxation+QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}
	
	if ( _THIS->growFarBounds < 1.1 )
	{
		_THIS->growFarBounds = 1.1;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->epsFlipping <= 0.0 )
	{
		_THIS->epsFlipping = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->numRegularisationSteps < 0 )
	{
		_THIS->numRegularisationSteps = 0;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->epsRegularisation < 0.0 )
	{
		_THIS->epsRegularisation = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->numRefinementSteps < 0 )
	{
		_THIS->numRefinementSteps = 0;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->epsIterRef < 0.0 )
	{
		_THIS->epsIterRef = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->epsLITests < 0.0 )
	{
		_THIS->epsLITests = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( _THIS->epsNZCTests < 0.0 )
	{
		_THIS->epsNZCTests = QPOASES_EPS;
		needToAdjust = BT_TRUE;
	}

	if ( needToAdjust == BT_TRUE )
		return THROWWARNING( RET_OPTIONS_ADJUSTED );

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue Options_print(	Options* _THIS
							)
{
	#ifndef __SUPPRESSANYOUTPUT__

	myStatic char myPrintfString[QPOASES_MAX_STRING_LENGTH];
	myStatic char info[QPOASES_MAX_STRING_LENGTH];

	qpOASES_myPrintf( "\n###################   qpOASES  --  QP OPTIONS   ##################\n" );
	qpOASES_myPrintf( "\n" );

	qpOASES_convertPrintLevelToString( _THIS->printLevel,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"printLevel                     =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n" );

	qpOASES_convertBooleanTypeToString( _THIS->enableRamping,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableRamping                  =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( _THIS->enableFarBounds,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableFarBounds                =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( _THIS->enableFlippingBounds,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableFlippingBounds           =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( _THIS->enableRegularisation,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableRegularisation           =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( _THIS->enableFullLITests,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableFullLITests              =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );
	
	qpOASES_convertBooleanTypeToString( _THIS->enableNZCTests,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableNZCTests                 =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableDriftCorrection          =  %d\n",_THIS->enableDriftCorrection );
	qpOASES_myPrintf( myPrintfString );
	
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableCholeskyRefactorisation  =  %d\n",_THIS->enableCholeskyRefactorisation );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertBooleanTypeToString( _THIS->enableEqualities,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"enableEqualities               =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n" );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"terminationTolerance           =  %e\n",_THIS->terminationTolerance );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"boundTolerance                 =  %e\n",_THIS->boundTolerance );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"boundRelaxation                =  %e\n",_THIS->boundRelaxation );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsNum                         =  %e\n",_THIS->epsNum );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsDen                         =  %e\n",_THIS->epsDen );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"maxPrimalJump                  =  %e\n",_THIS->maxPrimalJump );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"maxDualJump                    =  %e\n",_THIS->maxDualJump );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n" );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"initialRamping                 =  %e\n",_THIS->initialRamping );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"finalRamping                   =  %e\n",_THIS->finalRamping );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"initialFarBounds               =  %e\n",_THIS->initialFarBounds );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"growFarBounds                  =  %e\n",_THIS->growFarBounds );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_convertSubjectToStatusToString( _THIS->initialStatusBounds,info );
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"initialStatusBounds            =  %s\n",info );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsFlipping                    =  %e\n",_THIS->epsFlipping );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"numRegularisationSteps         =  %d\n",_THIS->numRegularisationSteps );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsRegularisation              =  %e\n",_THIS->epsRegularisation );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"numRefinementSteps             =  %d\n",_THIS->numRefinementSteps );
	qpOASES_myPrintf( myPrintfString );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsIterRef                     =  %e\n",_THIS->epsIterRef );
	qpOASES_myPrintf( myPrintfString );
	
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsLITests                     =  %e\n",_THIS->epsLITests );
	qpOASES_myPrintf( myPrintfString );
	
	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"epsNZCTests                    =  %e\n",_THIS->epsNZCTests );
	qpOASES_myPrintf( myPrintfString );

	qpOASES_myPrintf( "\n\n" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/



END_NAMESPACE_QPOASES


/*
 *	end of file
 */

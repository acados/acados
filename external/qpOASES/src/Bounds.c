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
 *	\file src/Bounds.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Implementation of the Bounds class designed to manage working sets of
 *	bounds within a QProblem.
 */


#include <qpOASES_e/Bounds.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	B o u n d s
 */
void BoundsCON( Bounds* _THIS, int _n )
{
	Bounds_init( _THIS,_n );
}


/*
 *	c o p y
 */
void BoundsCPY(	Bounds* FROM,
				Bounds* TO
				)
{
	int i;

	TO->n = FROM->n;
	TO->noLower = FROM->noLower;
	TO->noUpper = FROM->noUpper;

	if ( FROM->n != 0 )
	{
		for( i=0; i<TO->n; ++i )
		{
			TO->type[i]   = FROM->type[i];
			TO->status[i] = FROM->status[i];
		}
	}

	IndexlistCPY( &(FROM->freee),&(TO->freee) );
	IndexlistCPY( &(FROM->fixed),&(TO->fixed) );
}



/*
 *	i n i t
 */
returnValue Bounds_init(	Bounds* _THIS,
							int _n
							)
{
	int i;
	
	if ( _n < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	if ( _n >= 0 )
	{
		Indexlist_init( &(_THIS->freee),_n );
		Indexlist_init( &(_THIS->fixed),_n );
	}

	_THIS->n = _n;
	_THIS->noLower = BT_TRUE;
	_THIS->noUpper = BT_TRUE;

	assert( _THIS->n <= NVMAX );

	if ( _THIS->n > 0 )
	{
		for( i=0; i<_THIS->n; ++i )
		{
			_THIS->type[i]   = ST_UNKNOWN;
			_THIS->status[i] = ST_UNDEFINED;
		}
	}
	
	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p B o u n d
 */
returnValue Bounds_setupBound(	Bounds* _THIS, int number, SubjectToStatus _status
								)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= _THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Add bound index to respective index list. */
	switch ( _status )
	{
		case ST_INACTIVE:
			if ( Bounds_addIndex( _THIS,Bounds_getFree( _THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		case ST_LOWER:
			if ( Bounds_addIndex( _THIS,Bounds_getFixed( _THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		case ST_UPPER:
			if ( Bounds_addIndex( _THIS,Bounds_getFixed( _THIS ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		default:
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A l l F r e e
 */
returnValue Bounds_setupAllFree( Bounds* _THIS )
{
	return Bounds_setupAll( _THIS,ST_INACTIVE );
}


/*
 *	s e t u p A l l L o w e r
 */
returnValue Bounds_setupAllLower( Bounds* _THIS )
{
	return Bounds_setupAll( _THIS,ST_LOWER );
}


/*
 *	s e t u p A l l U p p e r
 */
returnValue Bounds_setupAllUpper( Bounds* _THIS )
{
	return Bounds_setupAll( _THIS,ST_UPPER );
}


/*
 *	m o v e F i x e d T o F r e e
 */
returnValue Bounds_moveFixedToFree( Bounds* _THIS, int number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= _THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of fixed variables to that of free ones. */
	if ( Bounds_removeIndex( _THIS,Bounds_getFixed( _THIS ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( Bounds_addIndex( _THIS,Bounds_getFree( _THIS ),number,ST_INACTIVE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	m o v e F r e e T o F i x e d
 */
returnValue Bounds_moveFreeToFixed(	Bounds* _THIS, int number, SubjectToStatus _status
									)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= _THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of free variables to that of fixed ones. */
	if ( Bounds_removeIndex( _THIS,Bounds_getFree( _THIS ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( Bounds_addIndex( _THIS,Bounds_getFixed( _THIS ),number,_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	f l i p F i x e d
 */
returnValue Bounds_flipFixed( Bounds* _THIS, int number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= _THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( _THIS->status != 0 )
		switch (_THIS->status[number])
		{
			case ST_LOWER: _THIS->status[number] = ST_UPPER; break;
			case ST_UPPER: _THIS->status[number] = ST_LOWER; break;
			default: return THROWERROR( RET_MOVING_BOUND_FAILED );
		}

	return SUCCESSFUL_RETURN;
}


/*
 *	s w a p F r e e
 */
returnValue Bounds_swapFree(	Bounds* _THIS, int number1, int number2
								)
{
	/* consistency check */
	if ( ( number1 < 0 ) || ( number1 >= _THIS->n ) || ( number2 < 0 ) || ( number2 >= _THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Swap index within indexlist of free variables. */
	return Bounds_swapIndex( _THIS,Bounds_getFree( _THIS ),number1,number2 );
}


/*
 *	s h i f t
 */
returnValue Bounds_shift(	Bounds* _THIS, int offset )
{
	int i;

	myStatic Indexlist shiftedFreee;
	myStatic Indexlist shiftedFixed;

	Indexlist_init( &shiftedFreee,_THIS->n );
	Indexlist_init( &shiftedFixed,_THIS->n );

	/* consistency check */
	if ( ( offset == 0 ) || ( _THIS->n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > _THIS->n/2 ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( ( _THIS->n % offset ) != 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Shift types and _THIS->status. */
	for( i=0; i<_THIS->n-offset; ++i )
	{
		Bounds_setType( _THIS,i,Bounds_getType( _THIS,i+offset ) );
		Bounds_setStatus( _THIS,i,Bounds_getStatus( _THIS,i+offset ) );
	}

	/* 2) Construct shifted index lists of free and fixed variables. */
	for( i=0; i<_THIS->n; ++i )
	{
		switch ( Bounds_getStatus( _THIS,i ) )
		{
			case ST_INACTIVE:
				if ( Indexlist_addNumber( &shiftedFreee,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_LOWER:
				if ( Indexlist_addNumber( &shiftedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_UPPER:
				if ( Indexlist_addNumber( &shiftedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			default:
				return THROWERROR( RET_SHIFTING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	IndexlistCPY( &shiftedFreee,&(_THIS->freee) );
	IndexlistCPY( &shiftedFixed,&(_THIS->fixed) );

	return SUCCESSFUL_RETURN;
}


/*
 *	r o t a t e
 */
returnValue Bounds_rotate( Bounds* _THIS, int offset )
{
	int i;
	myStatic SubjectToType   typeTmp[NVMAX];
	myStatic SubjectToStatus statusTmp[NVMAX];
	
	myStatic Indexlist rotatedFreee;
	myStatic Indexlist rotatedFixed;

	Indexlist_init( &rotatedFreee,_THIS->n );
	Indexlist_init( &rotatedFixed,_THIS->n );

	/* consistency check */
	if ( ( offset == 0 ) || ( offset == _THIS->n ) || ( _THIS->n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > _THIS->n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );


	/* 1) Rotate types and status. */
	for( i=0; i<offset; ++i )
	{
		typeTmp[i] = Bounds_getType( _THIS,i );
		statusTmp[i] = Bounds_getStatus( _THIS,i );
	}

	for( i=0; i<_THIS->n-offset; ++i )
	{
		Bounds_setType( _THIS,i,Bounds_getType( _THIS,i+offset ) );
		Bounds_setStatus( _THIS,i,Bounds_getStatus( _THIS,i+offset ) );
	}

	for( i=_THIS->n-offset; i<_THIS->n; ++i )
	{
		Bounds_setType( _THIS,i,typeTmp[i-_THIS->n+offset] );
		Bounds_setStatus( _THIS,i,statusTmp[i-_THIS->n+offset] );
	}

	/* 2) Construct shifted index lists of free and fixed variables. */
	for( i=0; i<_THIS->n; ++i )
	{
		switch ( Bounds_getStatus( _THIS,i ) )
		{
			case ST_INACTIVE:
				if ( Indexlist_addNumber( &rotatedFreee,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_LOWER:
				if ( Indexlist_addNumber( &rotatedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_UPPER:
				if ( Indexlist_addNumber( &rotatedFixed,i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			default:
				return THROWERROR( RET_ROTATING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	IndexlistCPY( &rotatedFreee,&(_THIS->freee) );
	IndexlistCPY( &rotatedFixed,&(_THIS->fixed) );

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue Bounds_print( Bounds* _THIS )
{
	#ifndef __SUPPRESSANYOUTPUT__

	myStatic char myPrintfString[QPOASES_MAX_STRING_LENGTH];

	int nFR = Bounds_getNFR( _THIS );
	int nFX = Bounds_getNFX( _THIS );

	int *FR_idx, *FX_idx;

	if ( _THIS->n == 0 )
		return SUCCESSFUL_RETURN;

	Indexlist_getNumberArray( Bounds_getFree( _THIS ),&FR_idx );
	Indexlist_getNumberArray( Bounds_getFixed( _THIS ),&FX_idx );

	snprintf( myPrintfString,QPOASES_MAX_STRING_LENGTH,"Bounds object comprising %d variables (%d free, %d fixed):\n",_THIS->n,nFR,nFX );
	qpOASES_myPrintf( myPrintfString );

	REFER_NAMESPACE_QPOASES qpOASES_printNI( FR_idx,nFR,"free " );
	REFER_NAMESPACE_QPOASES qpOASES_printNI( FX_idx,nFX,"fixed" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	s e t u p A l l
 */
returnValue Bounds_setupAll( Bounds* _THIS, SubjectToStatus _status )
{
	int i;

	/* 1) Place unbounded variables at the beginning of the index list of free variables. */
	for( i=0; i<_THIS->n; ++i )
	{
		if ( Bounds_getType( _THIS,i ) == ST_UNBOUNDED )
		{
			if ( Bounds_setupBound( _THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 2) Add remaining (i.e. bounded but possibly free) variables to the index list of free variables. */
	for( i=0; i<_THIS->n; ++i )
	{
		if ( Bounds_getType( _THIS,i ) == ST_BOUNDED )
		{
			if ( Bounds_setupBound( _THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 3) Place implicitly fixed variables at the end of the index list of free variables. */
	for( i=0; i<_THIS->n; ++i )
	{
		if ( Bounds_getType( _THIS,i ) == ST_EQUALITY )
		{
			if ( Bounds_setupBound( _THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 4) Moreover, add all bounds of unknown type. */
	for( i=0; i<_THIS->n; ++i )
	{
		if ( Bounds_getType( _THIS,i ) == ST_UNKNOWN || Bounds_getType( _THIS,i ) == ST_DISABLED )
		{
			if ( Bounds_setupBound( _THIS,i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	a d d I n d e x
 */
returnValue Bounds_addIndex(	Bounds* _THIS, Indexlist* const indexlist,
								int newnumber, SubjectToStatus newstatus
								)
{
	if ( _THIS->status != 0 )
	{
		/* consistency check */
		if ( _THIS->status[newnumber] == newstatus )
			return THROWERROR( RET_INDEX_ALREADY_OF_DESIRED_STATUS );

		_THIS->status[newnumber] = newstatus;
	}
	else
		return THROWERROR( RET_ADDINDEX_FAILED );

	if ( indexlist != 0 )
	{
		if ( Indexlist_addNumber( indexlist,newnumber ) == RET_INDEXLIST_EXCEEDS_MAX_LENGTH )
			return THROWERROR( RET_ADDINDEX_FAILED );
	}
	else
		return THROWERROR( RET_INVALID_ARGUMENTS );

	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e I n d e x
 */
returnValue Bounds_removeIndex(	Bounds* _THIS, Indexlist* const indexlist,
								int removenumber
								)
{
	if ( _THIS->status != 0 )
		_THIS->status[removenumber] = ST_UNDEFINED;
	else
		return THROWERROR( RET_REMOVEINDEX_FAILED );

	if ( indexlist != 0 )
	{
		if ( Indexlist_removeNumber( indexlist,removenumber ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_REMOVEINDEX_FAILED );
	}
	else
		return THROWERROR( RET_INVALID_ARGUMENTS );

	return SUCCESSFUL_RETURN;
}


/*
 *	s w a p I n d e x
 */
returnValue Bounds_swapIndex(	Bounds* _THIS, Indexlist* const indexlist,
								int number1, int number2
								)
{
	/* consistency checks */
	if ( _THIS->status != 0 )
	{
		if ( _THIS->status[number1] != _THIS->status[number2] )
			return THROWERROR( RET_SWAPINDEX_FAILED );
	}
	else
		return THROWERROR( RET_SWAPINDEX_FAILED );

	if ( number1 == number2 )
	{
		THROWWARNING( RET_NOTHING_TO_DO );
		return SUCCESSFUL_RETURN;
	}

	if ( indexlist != 0 )
	{
		if ( Indexlist_swapNumbers( indexlist,number1,number2 ) != SUCCESSFUL_RETURN )
			return THROWERROR( RET_SWAPINDEX_FAILED );
	}
	else
		return THROWERROR( RET_INVALID_ARGUMENTS );

	return SUCCESSFUL_RETURN;
}



END_NAMESPACE_QPOASES


/*
 *	end of file
 */

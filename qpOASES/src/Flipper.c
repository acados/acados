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
 *	\file src/Flipper.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Implementation of the Flipper class designed to manage working sets of
 *	constraints and bounds within a QProblem.
 */


#include <qpOASES_e/Flipper.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


 /*
 *	F l i p p e r
 */
void FlipperCON(	Flipper* _THIS,
					unsigned int _nV,
					unsigned int _nC
					)
{
	Flipper_init( _THIS,_nV,_nC );
}


/*
 *	c o p y
 */
void FlipperCPY(	Flipper* FROM,
					Flipper* TO
					)
{
	Flipper_set( TO, &(FROM->bounds),FROM->R,&(FROM->constraints),FROM->Q,FROM->T );
}


/*
 *	i n i t
 */
returnValue Flipper_init(	Flipper* _THIS,
							unsigned int _nV,
							unsigned int _nC
							)
{
	_THIS->nV = _nV;
	_THIS->nC = _nC;

	return SUCCESSFUL_RETURN;
}



/*
 *	g e t
 */
returnValue Flipper_get(	Flipper* _THIS,
							Bounds* const _bounds,
							real_t* const _R,
							Constraints* const _constraints,
							real_t* const _Q,
							real_t* const _T 
							)
{
	if ( _bounds != 0 )
		BoundsCPY( &(_THIS->bounds),_bounds );

	if ( _constraints != 0 )
		ConstraintsCPY( &(_THIS->constraints),_constraints );

	if ( _R != 0 )
		memcpy( _R,_THIS->R, (NVMAX*NVMAX)*sizeof(real_t) );

	if ( _Q != 0 )
		memcpy( _Q,_THIS->Q, (NVMAX*NVMAX)*sizeof(real_t) );

	if ( _T != 0 )
		memcpy( _T,_THIS->T, (NVCMIN*NVCMIN)*sizeof(real_t) );

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t
 */
returnValue Flipper_set(	Flipper* _THIS,
							const Bounds* const _bounds,
							const real_t* const _R,
							const Constraints* const _constraints,
							const real_t* const _Q,
							const real_t* const _T
							)
{
	if ( _bounds != 0 )
		BoundsCPY( (Bounds*)_bounds,&(_THIS->bounds) );

	if ( _constraints != 0 )
		ConstraintsCPY( (Constraints*)_constraints,&(_THIS->constraints) );

	if ( _R != 0 )
		memcpy( _THIS->R,_R, (NVMAX*NVMAX)*sizeof(real_t) );

	if ( _Q != 0 )
		memcpy( _THIS->Q,_Q, (NVMAX*NVMAX)*sizeof(real_t) );

	if ( _T != 0 )
		memcpy( _THIS->T,_T, (NVCMIN*NVCMIN)*sizeof(real_t) );

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/


unsigned int Flipper_getDimT( Flipper* _THIS )
{
	if ( _THIS->nV > _THIS->nC )
		return _THIS->nC * _THIS->nC;
	else
		return _THIS->nV * _THIS->nV;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */

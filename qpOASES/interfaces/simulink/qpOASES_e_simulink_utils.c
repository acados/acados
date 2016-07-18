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
 *	\file interfaces/simulink/qpOASES_e_simulink_utils.c
 *	\author Hans Joachim Ferreau
 *	\version 3.1
 *	\date 2007-2015
 *
 *	Collects utility functions for Interface to Simulink(R) that
 *	enables to call qpOASES as a C S function.
 *
 */


USING_NAMESPACE_QPOASES


/*
 *	i s N a N
 */
BooleanType isNaN( real_t val )
{
	if ( (( val <= 0.0 ) || ( val >= 0.0 )) == 0 )
		return BT_TRUE;
	else
		return BT_FALSE;
}


/*
 *	r e m o v e N a N s
 */
returnValue removeNaNs( real_t* const data, unsigned int dim )
{
	unsigned int i;

	if ( data == 0 )
		return RET_INVALID_ARGUMENTS;

	for ( i=0; i<dim; ++i )
		if ( isNaN(data[i]) == BT_TRUE )
			data[i] = QPOASES_INFTY;

	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e I n f s
 */
returnValue removeInfs( real_t* const data, unsigned int dim )
{
	unsigned int i;

	if ( data == 0 )
		return RET_INVALID_ARGUMENTS;

	for ( i=0; i<dim; ++i )
	{
		if ( data[i] < -QPOASES_INFTY )
			data[i] = -QPOASES_INFTY;

		if ( data[i] > QPOASES_INFTY )
			data[i] = QPOASES_INFTY;
	}

	return SUCCESSFUL_RETURN;
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
 *	end of file
 */

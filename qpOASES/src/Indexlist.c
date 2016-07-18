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
 *	\file src/Indexlist.c
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Implementation of the Indexlist class designed to manage index lists of
 *	constraints and bounds within a QProblem_SubjectTo.
 */


#include <qpOASES_e/Indexlist.h>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

/*
 *	I n d e x l i s t
 */
void IndexlistCON(	Indexlist* _THIS,
					int n
					)
{
	Indexlist_init( _THIS,n );
}


/*
 *	c o p y
 */
void IndexlistCPY(	Indexlist* FROM,
					Indexlist* TO
					)
{
	int i;
	
	if ( FROM != TO )
	{
		TO->length = FROM->length;
		TO->physicallength = FROM->physicallength;

		if ( FROM->number != 0 )
		{
			for( i=0; i<TO->physicallength; ++i )
				TO->number[i] = FROM->number[i];
			for( i=0; i<TO->physicallength; ++i )
				TO->iSort[i] = FROM->iSort[i];
		}
	}
}



/*
 *	i n i t
 */
returnValue Indexlist_init(	Indexlist* _THIS,
							int n
							)
{
	if ( n < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	_THIS->length = 0;
	_THIS->physicallength = n;

	assert( n <= NVCMAX );

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t N u m b e r A r r a y
 */
returnValue Indexlist_getNumberArray( Indexlist* _THIS, int** const numberarray )
{
	if (numberarray == 0)
		return THROWERROR( RET_INVALID_ARGUMENTS );

	*numberarray = _THIS->number;
	return SUCCESSFUL_RETURN;
}


/*
 *	g e t I S o r t A r r a y
 */
returnValue Indexlist_getISortArray( Indexlist* _THIS, int** const iSortArray )
{
	*iSortArray = _THIS->iSort;
	return SUCCESSFUL_RETURN;
}



/*
 *	g e t I n d e x
 */
int Indexlist_getIndex( Indexlist* _THIS, int givennumber )
{
	int myIndex = Indexlist_findInsert(_THIS,givennumber);
	return _THIS->number[_THIS->iSort[myIndex]] == givennumber ? _THIS->iSort[myIndex] : -1;
}


/*
 *	a d d N u m b e r
 */
returnValue Indexlist_addNumber( Indexlist* _THIS, int addnumber )
{
	int i, j;

	if ( _THIS->length >= _THIS->physicallength )
		return THROWERROR( RET_INDEXLIST_EXCEEDS_MAX_LENGTH );

	_THIS->number[_THIS->length] = addnumber;
	j = Indexlist_findInsert(_THIS,addnumber);
	for (i = _THIS->length; i > j+1; i--)
		_THIS->iSort[i] = _THIS->iSort[i-1];
	_THIS->iSort[j+1] = _THIS->length;
	++(_THIS->length);

	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e N u m b e r
 */
returnValue Indexlist_removeNumber( Indexlist* _THIS, int removenumber )
{
	int i;
	int idx = Indexlist_findInsert( _THIS,removenumber );
	int iSidx = _THIS->iSort[idx];

	/* nothing to be done if number is not contained in index set */
	if ( _THIS->number[iSidx] != removenumber )
		return SUCCESSFUL_RETURN;

	/* update sorted indices iSort first */
	for (i = 0; i < _THIS->length; i++)
		if (_THIS->iSort[i] > iSidx) _THIS->iSort[i]--;
	for (i = idx+1; i < _THIS->length; i++)
		_THIS->iSort[i-1] = _THIS->iSort[i];

	/* remove from numbers list */
	for( i=iSidx; i<_THIS->length-1; ++i )
		_THIS->number[i] = _THIS->number[i+1];
	_THIS->number[_THIS->length-1] = -1;

	--(_THIS->length);

	return SUCCESSFUL_RETURN;
}


/*
 *	s w a p N u m b e r s
 */
returnValue Indexlist_swapNumbers( Indexlist* _THIS, int number1, int number2 )
{
	int index1 = Indexlist_findInsert( _THIS,number1 );
	int index2 = Indexlist_findInsert( _THIS,number2 );
	int tmp;

	/* consistency check */
	if ( ( _THIS->number[_THIS->iSort[index1]] != number1 ) || ( _THIS->number[_THIS->iSort[index2]] != number2 ) )
		return THROWERROR( RET_INDEXLIST_CORRUPTED );

	/* swap numbers */
	tmp = _THIS->number[_THIS->iSort[index1]];
	_THIS->number[_THIS->iSort[index1]] = _THIS->number[_THIS->iSort[index2]];
	_THIS->number[_THIS->iSort[index2]] = tmp;
	/* swap sorting indices */
	tmp = _THIS->iSort[index1];
	_THIS->iSort[index1] = _THIS->iSort[index2];
	_THIS->iSort[index2] = tmp;

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/


int Indexlist_findInsert( Indexlist* _THIS, int i )
{
	int fst = 0, lst = _THIS->length-1, mid;

	/* quick check if index can be appended */
	if (_THIS->length == 0 || i < _THIS->number[_THIS->iSort[0]]) return -1;
	if (i >= _THIS->number[_THIS->iSort[_THIS->length-1]]) return _THIS->length-1;

	/* otherwise, perform bisection search */
	while (fst < lst - 1)
	{
		mid = (fst + lst) / 2;
		if (i >= _THIS->number[_THIS->iSort[mid]]) fst = mid;
		else lst = mid;
	}

	return fst;
}

END_NAMESPACE_QPOASES


/*
 *	end of file
 */

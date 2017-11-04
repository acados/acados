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
 *	\file examples/example1b.c
 *	\author Hans Joachim Ferreau
 *	\version 3.1embedded
 *	\date 2007-2015
 *
 *	Very simple example for testing qpOASES using the QProblemB class.
 */


#include <qpOASES_e.h>


/** Example for qpOASES main function using the QProblemB class. */
int main( )
{
	USING_NAMESPACE_QPOASES

	real_t stat, feas, cmpl;

	/* Setup data of first QP. */
	real_t H[2*2] = { 1.0, 0.0, 0.0, 0.5 };
	real_t g[2] = { 1.5, 1.0 };
	real_t lb[2] = { 0.5, -2.0 };
	real_t ub[2] = { 5.0, 2.0 };

	/* Setup data of second QP. */
	real_t g_new[2] = { 1.0, 1.5 };
	real_t lb_new[2] = { 0.0, -1.0 };
	real_t ub_new[2] = { 5.0, -0.5 };


	/* Setting up QProblemB object. */
	QProblemB *example = QProblemB_createMemory(2);
	static Options options;

	int nWSR = 10;
	real_t xOpt[2];
	real_t yOpt[2];

	QProblemBCON( example,2,HST_UNKNOWN );
	Options_setToDefault( &options );
	/* options.enableFlippingBounds = BT_FALSE; */
	options.initialStatusBounds = ST_INACTIVE;
	options.numRefinementSteps = 1;
	/* options.enableCholeskyRefactorisation = 1; */
	QProblemB_setOptions( example,options );


	/* Solve first QP. */
	nWSR = 10;
	QProblemB_init( example,H,g,lb,ub, &nWSR,0 );

	/* Get and print solution of second QP. */
	QProblemB_getPrimalSolution( example,xOpt );
	QProblemB_getDualSolution( example,yOpt );
	printf( "\nxOpt = [ %e, %e ];  objVal = %e\n\n", xOpt[0],xOpt[1],QProblemB_getObjVal(example) );

	qpOASES_getKktViolation( 2,1, H,g,NULL,lb,ub,NULL,NULL, xOpt,yOpt, &stat,&feas,&cmpl );
	printf("KKT violations:\n\n");
	printf("stat = %e, feas = %e, cmpl = %e\n\n", stat, feas, cmpl);


	/* Solve second QP. */
	nWSR = 10;
	QProblemB_hotstart( example,g_new,lb_new,ub_new, &nWSR,0 );

	/* Get and print solution of second QP. */
	QProblemB_getPrimalSolution( example,xOpt );
	QProblemB_getDualSolution( example,yOpt );
	printf( "\nxOpt = [ %e, %e ];  objVal = %e\n\n", xOpt[0],xOpt[1],QProblemB_getObjVal(example) );

	qpOASES_getKktViolation( 2,1, H,g_new,NULL,lb_new,ub_new,NULL,NULL, xOpt,yOpt, &stat,&feas,&cmpl );
	printf("KKT violations:\n\n");
	printf("stat = %e, feas = %e, cmpl = %e\n\n", stat, feas, cmpl);

	free(example);

	return 0;
}


/*
 *	end of file
 */

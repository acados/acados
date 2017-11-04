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
 *	\file examples/exampleLP.c
 *	\author Hans Joachim Ferreau
 *	\version 3.1embedded
 *	\date 2008-2015
 *
 *	Very simple example for solving a LP sequence using qpOASES.
 */


#include <qpOASES_e.h>


/** Example for qpOASES main function solving LPs. */
int main( )
{
	USING_NAMESPACE_QPOASES

	real_t stat, feas, cmpl;

	/* Setup data of first LP. */
	real_t A[1*2] = { 1.0, 1.0 };
	real_t g[2] = { 1.5, 1.0 };
	real_t lb[2] = { 0.5, -2.0 };
	real_t ub[2] = { 5.0, 2.0 };
	real_t lbA[1] = { -1.0 };
	real_t ubA[1] = { 2.0 };

	/* Setup data of second LP. */
	real_t g_new[2] = { 1.0, 1.5 };
	real_t lb_new[2] = { 0.0, -1.0 };
	real_t ub_new[2] = { 5.0, -0.5 };
	real_t lbA_new[1] = { -2.0 };
	real_t ubA_new[1] = { 1.0 };


	/* Setting up QProblem object with zero Hessian matrix. */
	QProblem *example = QProblem_createMemory(2, 1);
	static Options options;

	int nWSR;
	real_t xOpt[2];
	real_t yOpt[2+1];

	QProblemCON( example,2,1,HST_ZERO );
	Options_setToDefault( &options );
	QProblem_setOptions( example,options );

	/* Solve first LP. */
	nWSR = 10;
	QProblem_init( example, 0,g,A,lb,ub,lbA,ubA, &nWSR,0 );

	/* Solve second LP. */
	nWSR = 10;
	QProblem_hotstart( example,g_new,lb_new,ub_new,lbA_new,ubA_new, &nWSR,0 );

	/* Get and print solution of second LP. */
	QProblem_getPrimalSolution( example,xOpt );
	QProblem_getDualSolution( example,yOpt );
	printf( "\nxOpt = [ %e, %e ];  objVal = %e\n\n", xOpt[0],xOpt[1],QProblem_getObjVal(example) );

	qpOASES_getKktViolation( 2,1, NULL,g_new,A,lb_new,ub_new,lbA_new,ubA_new, xOpt,yOpt, &stat,&feas,&cmpl );
	printf("KKT violations:\n\n");
	printf("stat = %e, feas = %e, cmpl = %e\n\n", stat, feas, cmpl);

	free(example);

	return 0;
}


/*
 *	end of file
 */

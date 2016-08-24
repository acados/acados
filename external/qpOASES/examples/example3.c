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
 *	\file examples/example3.c
 *	\author Hans Joachim Ferreau
 *	\version 3.1embedded
 *	\date 2008-2015
 *
 *	Example demonstrating usage of qpOASES for solving a QP sequence of the
 *	Online QP Benchmark Collection. In order to run it, you have to download
 *	"Example 02" from from http://www.qpOASES.org/onlineQP/ and store it into 
 *	the directory bin/chain80w/.
 */



#include <qpOASES_e.h>


/** Example for qpOASES main function using the OQP interface. */
int main( )
{
	USING_NAMESPACE_QPOASES

	/* 1) Define benchmark arguments. */
	BooleanType isSparse = BT_FALSE;
	BooleanType useHotstarts = BT_TRUE;
	
	int maxAllowedNWSR = 600;
	real_t maxNWSR, avgNWSR, maxCPUtime, avgCPUtime;
	real_t maxStationarity, maxFeasibility, maxComplementarity;

	static Options options;
	Options_setToMPC( &options );
	options.printLevel = PL_LOW;

	/* 2) Run benchmark. */
	if ( runOQPbenchmark(	"./chain80w/",
							isSparse,useHotstarts,
							&options,maxAllowedNWSR,
							&maxNWSR,&avgNWSR,&maxCPUtime,&avgCPUtime,
							&maxStationarity,&maxFeasibility,&maxComplementarity
							) != SUCCESSFUL_RETURN )
	{
		qpOASES_myPrintf( "In order to run this example, you need to download example no. 02\nfrom the Online QP Benchmark Collection website first!\n" );
		fprintf( stderr,"error\n" );
		return -1;
	}

	/* 3) Print results. */
	fprintf( stderr,"\n\n" );
	fprintf( stderr,"OQP Benchmark Results:\n" );
	fprintf( stderr,"======================\n\n" );
	fprintf( stderr,"maximum stationary error:     %.3e\n",maxStationarity );
	fprintf( stderr,"maximum feasibility error:    %.3e\n",maxFeasibility );
	fprintf( stderr,"maximum complementary error:  %.3e\n",maxComplementarity );
	fprintf( stderr,"\n" );
	fprintf( stderr,"maximum CPU time:             %.3f milliseconds\n\n",1000.0*maxCPUtime );

	return 0;
}


/*
 *	end of file
 */

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

#define _SVID_SOURCE


#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <qpOASES_e.h>
#include <qpOASES_e/UnitTesting.h>


/** Try to solve a list of or all OQP examples in testing/problems */
int main(int argc, char *argv[])
{
	const real_t TOL = 1e-5;

	/* 1) Define benchmark arguments. */
	BooleanType isSparse = BT_FALSE;
	BooleanType useHotstarts = BT_FALSE;
	static Options options;

	int maxAllowedNWSR;
	real_t maxNWSR, avgNWSR, maxCPUtime, avgCPUtime;
	real_t maxStationarity, maxFeasibility, maxComplementarity;

	int scannedDir = 0;
	int nfail = 0, npass = 0;
	int nproblems, i;
	struct dirent **namelist;
	char resstr[200], OQPproblem[200];
	char *problem;
	returnValue returnvalue;

	Options_setToDefault( &options );
	/*options.enableFlippingBounds = BT_FALSE;*/
	/*Options_setToReliable( &options );*/
	Options_setToMPC( &options );
	/*options.printLevel = PL_DEBUG_ITER;*/
	options.printLevel = PL_LOW;
	/*options.enableRamping = BT_FALSE;*/
	/*options.enableFarBounds = BT_FALSE;*/
	
	if (argc == 1)
	{
		/* 2a) Scan problem directory */
		nproblems = scandir("../testing/c/data/problems", &namelist, NULL, alphasort);
		if (nproblems <= 0)
		{
			qpOASES_myPrintf( "No test problems found!\n" );
			return -1;
		}
		scannedDir = 1;
	}
	else
	{
		/* 2b) Use problem list given by arguments */
		nproblems = argc - 1;
		scannedDir = 0;
	}

	/* 3) Run benchmark. */
	printf("%10s %9s %9s %9s %6s  %-12s\n", "problem", "stat",
			"feas", "compl", "nWSR", "result");
	for (i = 0; i < nproblems; i++)
	{
		if (scannedDir)
		{
			/* skip special directories and zip file cuter.*bz2 */
			if (namelist[i]->d_name[0] == '.' || namelist[i]->d_name[0] == 'c')
			{
				free(namelist[i]);
				continue;
			}
			problem = namelist[i]->d_name;
		}
		else
		{
			problem = argv[i+1];
		}

		fprintf(stdout, "%-10s ", problem);
		fflush(stdout);

		snprintf(OQPproblem, 199, "../testing/c/data/problems/%s/", problem);
		maxCPUtime = 300.0;
		maxAllowedNWSR = 3500;
		returnvalue = runOQPbenchmark(	OQPproblem,
										isSparse,useHotstarts,
										&options,maxAllowedNWSR,
										&maxNWSR,&avgNWSR,&maxCPUtime,&avgCPUtime,
										&maxStationarity,&maxFeasibility,&maxComplementarity
										);
		if (returnvalue	== SUCCESSFUL_RETURN
				&& maxStationarity < TOL
				&& maxFeasibility < TOL
				&& maxComplementarity < TOL)
		{
			npass++;
			strncpy(resstr, "pass", 199);
		}
		else
		{
			nfail++;
			snprintf (resstr, 199, "fail (%d)", returnvalue);
		}
		fprintf(stdout, "%9.2e %9.2e %9.2e %6d  %-12s\n", maxStationarity,
				maxFeasibility, maxComplementarity, (int)maxNWSR, resstr);

		if (scannedDir) free(namelist[i]);
	}
	if (scannedDir) free(namelist);

	/* 4) Print results. */
	printf("\n\n" );
	printf("Testbench results:\n" );
	printf("======================\n\n" );
	printf("Pass:  %3d\n", npass);
	printf("Fail:  %3d\n", nfail);
	printf("Ratio: %5.1f%%\n", 100.0 * (real_t)npass / (real_t)(npass+nfail));
	printf("\n" );

	return 0;
}


/*
 *	end of file
 */

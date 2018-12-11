// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados/sim/sim_common.h"
#include "acados_c/sim_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	long long *ptr;
	int ii;



	/* RHS */

	// model

	int nu = mxGetScalar( mxGetField( prhs[0], 0, "nu" ) );
//	mexPrintf("\n%d\n", nu);
	int nx = mxGetScalar( mxGetField( prhs[0], 0, "nx" ) );
//	mexPrintf("\n%d\n", nx);


	// opts_struct

	int num_stages = mxGetScalar( mxGetField( prhs[1], 0, "num_stages" ) );
	int num_steps = mxGetScalar( mxGetField( prhs[1], 0, "num_steps" ) );
	bool sens_forw = mxGetScalar( mxGetField( prhs[1], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
	char *scheme = mxArrayToString( mxGetField( prhs[1], 0, "scheme" ) );
//	mexPrintf("\n%s\n", scheme);



	/* LHS */

	// field names of output struct
	char *fieldnames[6];
	fieldnames[0] = (char*)mxMalloc(50);
	fieldnames[1] = (char*)mxMalloc(50);
	fieldnames[2] = (char*)mxMalloc(50);
	fieldnames[3] = (char*)mxMalloc(50);
	fieldnames[4] = (char*)mxMalloc(50);
	fieldnames[5] = (char*)mxMalloc(50);

	memcpy(fieldnames[0],"config",sizeof("config"));
	memcpy(fieldnames[1],"dims",sizeof("dims"));
	memcpy(fieldnames[2],"opts",sizeof("opts"));
	memcpy(fieldnames[3],"in",sizeof("in"));
	memcpy(fieldnames[4],"out",sizeof("out"));
	memcpy(fieldnames[5],"solver",sizeof("solver"));

	// create output struct
	plhs[0] = mxCreateStructMatrix(1, 1, 6, (const char **) fieldnames);

	mxFree( fieldnames[0] );
	mxFree( fieldnames[1] );
	mxFree( fieldnames[2] );
	mxFree( fieldnames[3] );
	mxFree( fieldnames[4] );
	mxFree( fieldnames[5] );



	/* plan & config */
	sim_solver_plan plan;

	if(!strcmp(scheme, "erk"))
		{
//		mexPrintf("\n%s\n", scheme);
		plan.sim_solver = ERK;
		}
	else if(!strcmp(scheme, "irk"))
		{
//		mexPrintf("\n%s\n", scheme);
		plan.sim_solver = IRK;
		}
	else
		{
		mexPrintf("\nscheme not supported %s\n", scheme);
		return;
		}

	sim_solver_config *config = sim_config_create(plan);


	/* dims */
	void *dims = sim_dims_create(config);
	sim_dims_set(config, dims, "nx", &nx);
	sim_dims_set(config, dims, "nu", &nu);


	/* opts */
	sim_rk_opts *opts = sim_opts_create(config, dims);
	sim_rk_opts_set(opts, "num_stages", &num_stages);
	sim_rk_opts_set(opts, "num_steps", &num_steps);
	sim_rk_opts_set(opts, "sens_forw", &sens_forw);

	/* in */
	sim_in *in = sim_in_create(config, dims);
	if(sens_forw==true)
		{
//		mexPrintf("\nsens forw true!\n");
		double *Sx = calloc(nx*nx, sizeof(double));
		for(ii=0; ii<nx; ii++)
			Sx[ii*(nx+1)] = 1.0;
		double *Su = calloc(nx*nu, sizeof(double));
//		d_print_mat(nx, nx, Sx, nx);
//		d_print_mat(nx, nu, Su, nx);
		sim_in_set(config, dims, in, "Sx", Sx);
		sim_in_set(config, dims, in, "Su", Su);
		free(Sx);
		free(Su);
		}


	/* out */
	sim_out *out = sim_out_create(config, dims);


	/* solver */
	sim_solver *solver = sim_solver_create(config, dims, opts);



	/* populate output struct */

	// config
	mxArray *config_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	ptr = mxGetData(config_mat);
	ptr[0] = (long long) config;
	mxSetField(plhs[0], 0, "config", config_mat);

	// dims
	mxArray *dims_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	ptr = mxGetData(dims_mat);
	ptr[0] = (long long) dims;
	mxSetField(plhs[0], 0, "dims", dims_mat);

	// opts
	mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	ptr = mxGetData(opts_mat);
	ptr[0] = (long long) opts;
	mxSetField(plhs[0], 0, "opts", opts_mat);

	// in
	mxArray *in_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	ptr = mxGetData(in_mat);
	ptr[0] = (long long) in;
	mxSetField(plhs[0], 0, "in", in_mat);

	// out
	mxArray *out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	ptr = mxGetData(out_mat);
	ptr[0] = (long long) out;
	mxSetField(plhs[0], 0, "out", out_mat);

	// solver
	mxArray *solver_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	ptr = mxGetData(solver_mat);
	ptr[0] = (long long) solver;
	mxSetField(plhs[0], 0, "solver", solver_mat);



	/* return */
	return;

	}

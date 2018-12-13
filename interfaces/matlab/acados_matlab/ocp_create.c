// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
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

	bool
		set_N = false,
		set_Ts = false,
		set_nx = false,
		set_nu = false,
		set_nbx = false,
		set_nbu = false,

	int N, nu, nx, nbx, nbu;
	double Ts;

	if(mxGetField( prhs[0], 0, "N" )!=NULL)
		{
		set_N = true;
		N = mxGetScalar( mxGetField( prhs[0], 0, "N" ) );
		}
	if(mxGetField( prhs[0], 0, "Ts" )!=NULL)
		{
		set_Ts = true;
		Ts = mxGetScalar( mxGetField( prhs[0], 0, "Ts" ) );
		}
	if(mxGetField( prhs[0], 0, "nx" )!=NULL)
		{
		set_nx = true;
		nx = mxGetScalar( mxGetField( prhs[0], 0, "nx" ) );
		}
	if(mxGetField( prhs[0], 0, "nu" )!=NULL)
		{
		set_nu = true;
		nu = mxGetScalar( mxGetField( prhs[0], 0, "nu" ) );
		}
	if(mxGetField( prhs[0], 0, "nbx" )!=NULL)
		{
		set_nbx = true;
		nbx = mxGetScalar( mxGetField( prhs[0], 0, "nbx" ) );
		}
	if(mxGetField( prhs[0], 0, "nbu" )!=NULL)
		{
		set_nbu = true;
		nbu = mxGetScalar( mxGetField( prhs[0], 0, "nbu" ) );
		}


	// opts_struct

//	int num_stages = mxGetScalar( mxGetField( prhs[1], 0, "num_stages" ) );
//	int num_steps = mxGetScalar( mxGetField( prhs[1], 0, "num_steps" ) );
//	bool sens_forw = mxGetScalar( mxGetField( prhs[1], 0, "sens_forw" ) );
//	mexPrintf("\n%d\n", sens_forw);
//	char *scheme = mxArrayToString( mxGetField( prhs[1], 0, "scheme" ) );
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
	ocp_nlp_plan *plan = ocp_nlp_plan_create(N);

	plan->nlp_solver = SQP;
//	plan->nlp_solver = SQP_RTI;
	
	ocp_nlp_plan_destroy(plan);
	
#if 0
	for(ii=0; ii<

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

	sim_config *config = sim_config_create(plan);


	/* dims */
	void *dims = sim_dims_create(config);
	if(set_nx)
		sim_dims_set(config, dims, "nx", &nx);
	if(set_nu)
		sim_dims_set(config, dims, "nu", &nu);


	/* opts */
	sim_opts *opts = sim_opts_create(config, dims);
	sim_opts_set(opts, "num_stages", &num_stages);
	sim_opts_set(opts, "num_steps", &num_steps);
	sim_opts_set(opts, "sens_forw", &sens_forw);

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
	if(set_T)
		sim_in_set(config, dims, in, "T", &T);
	if(set_x)
		sim_in_set(config, dims, in, "x", x);
	if(set_u)
		sim_in_set(config, dims, in, "u", u);


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
#endif



	/* return */
	return;

	}


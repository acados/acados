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
	long long *l_ptr;

	int ii;



	/* RHS */

	// model

	bool
		set_N = false,
		set_Ts = false,
		set_nx = false,
		set_nu = false,
		set_nbx = false,
		set_nbu = false;

	int N, nu, nx, nbx, nbu;
	double Ts;

	if(mxGetField( prhs[0], 0, "N" )!=NULL)
		{
		set_N = true;
		N = mxGetScalar( mxGetField( prhs[0], 0, "N" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: N not set!\n");
		return;
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
	else
		{
		mexPrintf("\nerror: ocp_create: nx not set!\n");
		return;
		}
	if(mxGetField( prhs[0], 0, "nu" )!=NULL)
		{
		set_nu = true;
		nu = mxGetScalar( mxGetField( prhs[0], 0, "nu" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: nu not set!\n");
		return;
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
	char *sim_scheme = mxArrayToString( mxGetField( prhs[1], 0, "sim_scheme" ) );
//	mexPrintf("\n%s\n", scheme);



	/* LHS */

	// field names of output struct
	char *fieldnames[6];
	fieldnames[0] = (char*) mxMalloc(50);
	fieldnames[1] = (char*) mxMalloc(50);
	fieldnames[2] = (char*) mxMalloc(50);
	fieldnames[3] = (char*) mxMalloc(50);
	fieldnames[4] = (char*) mxMalloc(50);
	fieldnames[5] = (char*) mxMalloc(50);

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

	// nlp solver
	plan->nlp_solver = SQP;
//	plan->nlp_solver = SQP_RTI;
	
	// cost
	for(ii=0; ii<=N; ii++)
		{
		plan->nlp_cost[ii] = LINEAR_LS;
//		plan->nlp_cost[ii] = NONLINEAR_LS;
		}
	
	// dynamics
	for(ii=0; ii<N; ii++)
		{
//		plan->nlp_dynamics[ii] = DISCRETE_MODEL;
		plan->nlp_dynamics[ii] = CONTINUOUS_MODEL;
		// TODO opts for that!!!
		if(!strcmp(sim_scheme, "erk"))
			{
			plan->sim_solver_plan[ii].sim_solver = ERK;
			}
		else if(!strcmp(sim_scheme, "irk"))
			{
			plan->sim_solver_plan[ii].sim_solver = IRK;
			}
		else
			{
			mexPrintf("\nsim_scheme not supported %s\n", sim_scheme);
			return;
			}
//		plan->sim_solver_plan[ii].sim_solver = LIFTED_IRK;
//		plan->sim_solver_plan[ii].sim_solver = GNSF;
		}
	
	for(ii=0; ii<=N; ii++)
		{
		plan->nlp_constraints[ii] = BGH;
//		plan->nlp_constraints[ii] = BGHP;
		}

    ocp_nlp_config *config = ocp_nlp_config_create(*plan);


	ocp_nlp_plan_destroy(plan);
	


	/* dims */

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	// allocate tmp
	int *i_ptr = (int *) malloc((N+1)*sizeof(int));
	// nx
	for(ii=0; ii<=N; ii++)
		i_ptr[ii] = nx;
	ocp_nlp_dims_set_opt_vars(config, dims, "nx", i_ptr);
	// nu
	for(ii=0; ii<N; ii++)
		i_ptr[ii] = nu;
	i_ptr[N] = 0;
	ocp_nlp_dims_set_opt_vars(config, dims, "nu", i_ptr);
	// free tmp
	free(i_ptr);
	// nbx
	ocp_nlp_dims_set_constraints(config, dims, 0, "nbx", &nx);
	if(set_nbx)
		{
		for(ii=1; ii<=N; ii++)
			ocp_nlp_dims_set_constraints(config, dims, ii, "nbx", &nbx);
		}
	// nbu
	if(set_nbu)
		{
		for(ii=0; ii<N; ii++)
			ocp_nlp_dims_set_constraints(config, dims, ii, "nbu", &nbu);
		}
			


	/* opts */
//	void *opts = ocp_nlp_opts_create(config, dims);
//	ocp_nlp_opts_destroy(opts);
#if 0
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
#endif



	/* populate output struct */

	// config
	mxArray *config_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(config_mat);
	l_ptr[0] = (long long) config;
	mxSetField(plhs[0], 0, "config", config_mat);

	// dims
	mxArray *dims_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(dims_mat);
	l_ptr[0] = (long long) dims;
	mxSetField(plhs[0], 0, "dims", dims_mat);

	// opts
//	mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
//	l_ptr = mxGetData(opts_mat);
//	l_ptr[0] = (long long) opts;
//	mxSetField(plhs[0], 0, "opts", opts_mat);

#if 0
	// in
	mxArray *in_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(in_mat);
	l_ptr[0] = (long long) in;
	mxSetField(plhs[0], 0, "in", in_mat);

	// out
	mxArray *out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(out_mat);
	l_ptr[0] = (long long) out;
	mxSetField(plhs[0], 0, "out", out_mat);

	// solver
	mxArray *solver_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(solver_mat);
	l_ptr[0] = (long long) solver;
	mxSetField(plhs[0], 0, "solver", solver_mat);
#endif



	/* return */
	return;

	}


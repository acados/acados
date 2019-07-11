// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados/sim/sim_common.h" // remove ???
#include "acados_c/sim_interface.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin sim_create\n");

	// sizeof(long long) == sizeof(void *) = 64 !!!
	int ii;

	long long *l_ptr;
	char *c_ptr;



	/* RHS */

	// opts_struct

	//
	int num_stages = mxGetScalar( mxGetField( prhs[1], 0, "num_stages" ) );
	//
	int num_steps = mxGetScalar( mxGetField( prhs[1], 0, "num_steps" ) );
	//
	bool sens_forw = false;
	c_ptr = mxArrayToString( mxGetField( prhs[1], 0, "sens_forw" ) );
	if (!strcmp(c_ptr, "true"))
		{
		sens_forw = true;
		}
//	mexPrintf("\n%d\n", sens_forw);
	//
	bool sens_adj = false;
	c_ptr = mxArrayToString( mxGetField( prhs[1], 0, "sens_adj" ) );
	if (!strcmp(c_ptr, "true"))
		{
		sens_adj = true;
		}
//	mexPrintf("\n%d\n", sens_adj);
	bool sens_hess = false;
	c_ptr = mxArrayToString( mxGetField( prhs[1], 0, "sens_hess" ) );
	if (!strcmp(c_ptr, "true"))
		{
		sens_hess = true;
		}
//	mexPrintf("\n%d\n", sens_hess);
	//
	char *method = mxArrayToString( mxGetField( prhs[1], 0, "method" ) );
//	mexPrintf("\n%s\n", method);



	// model

	int nu;				bool set_nu = false;
	int nx;				bool set_nx = false;
	double T;			bool set_T = false;
	double *x;			bool set_x = false;
	double *u;			bool set_u = false;
	double *seed_adj;	bool set_seed_adj = false;
	// gnsf stuff
	int gnsf_nx1;		bool set_gnsf_nx1 = false;
	int gnsf_nz1;		bool set_gnsf_nz1 = false;
	int gnsf_nuhat;		bool set_gnsf_nuhat = false;
	int gnsf_ny;		bool set_gnsf_ny = false;
	int gnsf_nout;		bool set_gnsf_nout = false;

	if(mxGetField( prhs[0], 0, "dim_nx" )!=NULL)
		{
		set_nx = true;
		nx = mxGetScalar( mxGetField( prhs[0], 0, "dim_nx" ) );
		}
	if(mxGetField( prhs[0], 0, "dim_nu" )!=NULL)
		{
		set_nu = true;
		nu = mxGetScalar( mxGetField( prhs[0], 0, "dim_nu" ) );
		}
	if(mxGetField( prhs[0], 0, "T" )!=NULL)
		{
		set_T = true;
		T = mxGetScalar( mxGetField( prhs[0], 0, "T" ) );
		}
	if(mxGetField( prhs[0], 0, "x" )!=NULL)
		{
		set_x = true;
		x = mxGetPr( mxGetField( prhs[0], 0, "x" ) );
		}
	if(mxGetField( prhs[0], 0, "u" )!=NULL)
		{
		set_u = true;
		u = mxGetPr( mxGetField( prhs[0], 0, "u" ) );
		}
	if(mxGetField( prhs[0], 0, "seed_adj" )!=NULL)
		{
		set_seed_adj = true;
		seed_adj = mxGetPr( mxGetField( prhs[0], 0, "seed_adj" ) );
		}
	// gnsf stuff
	if(!strcmp(method, "irk_gnsf"))
		{
		if(mxGetField( prhs[0], 0, "dim_gnsf_nx1" )!=NULL)
			{
			set_gnsf_nx1 = true;
			gnsf_nx1 = mxGetScalar( mxGetField( prhs[0], 0, "dim_gnsf_nx1" ) );
			}
		if(mxGetField( prhs[0], 0, "dim_gnsf_nz1" )!=NULL)
			{
			set_gnsf_nz1 = true;
			gnsf_nz1 = mxGetScalar( mxGetField( prhs[0], 0, "dim_gnsf_nz1" ) );
			}
		if(mxGetField( prhs[0], 0, "dim_gnsf_nuhat" )!=NULL)
			{
			set_gnsf_nuhat = true;
			gnsf_nuhat = mxGetScalar( mxGetField( prhs[0], 0, "dim_gnsf_nuhat" ) );
			}
		if(mxGetField( prhs[0], 0, "dim_gnsf_ny" )!=NULL)
			{
			set_gnsf_ny = true;
			gnsf_ny = mxGetScalar( mxGetField( prhs[0], 0, "dim_gnsf_ny" ) );
			}
		if(mxGetField( prhs[0], 0, "dim_gnsf_nout" )!=NULL)
			{
			set_gnsf_nout = true;
			gnsf_nout = mxGetScalar( mxGetField( prhs[0], 0, "dim_gnsf_nout" ) );
			}
		}


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

	if(!strcmp(method, "erk"))
		{
		plan.sim_solver = ERK;
		}
	else if(!strcmp(method, "irk"))
		{
		plan.sim_solver = IRK;
		}
	else if(!strcmp(method, "irk_gnsf"))
		{
		plan.sim_solver = GNSF;
		}
	else
		{
		mexPrintf("\nsim_create: method not supported %s\n", method);
		return;
		}

	sim_config *config = sim_config_create(plan);


	/* dims */
	void *dims = sim_dims_create(config);
	if(set_nx)
		sim_dims_set(config, dims, "nx", &nx);
	if(set_nu)
		sim_dims_set(config, dims, "nu", &nu);
	if(!strcmp(method, "irk_gnsf"))
		{
		if(set_gnsf_nx1)
			sim_dims_set(config, dims, "gnsf_nx1", &gnsf_nx1);
		if(set_gnsf_nz1)
			sim_dims_set(config, dims, "gnsf_nz1", &gnsf_nz1);
		if(set_gnsf_nuhat)
			sim_dims_set(config, dims, "gnsf_nuhat", &gnsf_nuhat);
		if(set_gnsf_ny)
			sim_dims_set(config, dims, "gnsf_ny", &gnsf_ny);
		if(set_gnsf_nout)
			sim_dims_set(config, dims, "gnsf_nout", &gnsf_nout);
		}


	/* opts */
	sim_opts *opts = sim_opts_create(config, dims);
	sim_opts_set(config, opts, "num_stages", &num_stages);
	sim_opts_set(config, opts, "num_steps", &num_steps);
	sim_opts_set(config, opts, "sens_forw", &sens_forw);
	sim_opts_set(config, opts, "sens_adj", &sens_adj);
	sim_opts_set(config, opts, "sens_hess", &sens_hess);


	/* in */
	sim_in *in = sim_in_create(config, dims);
	if(sens_forw==true) // | sens_hess==true ???
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
		{
		sim_in_set(config, dims, in, "T", &T);
		}
	if(set_x)
		{
		sim_in_set(config, dims, in, "x", x);
		}
	if(set_u)
		{
		sim_in_set(config, dims, in, "u", u);
		}
	if(set_seed_adj)
		{
		sim_in_set(config, dims, in, "S_adj", seed_adj);
		}


	/* out */
	sim_out *out = sim_out_create(config, dims);


	/* solver */
	sim_solver *solver = sim_solver_create(config, dims, opts);



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
	mxArray *opts_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(opts_mat);
	l_ptr[0] = (long long) opts;
	mxSetField(plhs[0], 0, "opts", opts_mat);

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



	/* return */
	return;

	}

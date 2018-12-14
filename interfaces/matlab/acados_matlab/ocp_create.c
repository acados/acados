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
	int *i_ptr;

	int ii, jj, idx;



	/* RHS */

	// model

	// dims
	double T;		bool set_T = false;
	int nx;			bool set_nx = false;
	int nu;			bool set_nu = false;
	int nyl;		bool set_nyl = false;
	int nym;		bool set_nym = false;
	int nbx;		bool set_nbx = false;
	int nbu;		bool set_nbu = false;
	int ng;			bool set_ng = false;
	// cost
	double *Vul;	bool set_Vul = false;
	double *Vxl;	bool set_Vxl = false;
	double *Vxm;	bool set_Vxm = false;
	double *Wl;		bool set_Wl = false;
	double *Wm;		bool set_Wm = false;
	double *yrl;	bool set_yrl = false;
	double *yrm;	bool set_yrm = false;
	// constraints
	double *x0;		bool set_x0 = false;
	double *Jbx;	bool set_Jbx = false;
	double *lbx;	bool set_lbx = false;
	double *ubx;	bool set_ubx = false;
	double *Jbu;	bool set_Jbu = false;
	double *lbu;	bool set_lbu = false;
	double *ubu;	bool set_ubu = false;
	double *C;		bool set_C = false;
	double *D;		bool set_D = false;
	double *lg;		bool set_lg = false;
	double *ug;		bool set_ug = false;

	// T
	if(mxGetField( prhs[0], 0, "T" )!=NULL)
		{
		set_T = true;
		T = mxGetScalar( mxGetField( prhs[0], 0, "T" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: T not set!\n");
		return;
		}
	// nx
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
	// nu
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
	// nyl
	if(mxGetField( prhs[0], 0, "nyl" )!=NULL)
		{
		set_nyl = true;
		nyl = mxGetScalar( mxGetField( prhs[0], 0, "nyl" ) );
		}
	// nym
	if(mxGetField( prhs[0], 0, "nym" )!=NULL)
		{
		set_nym = true;
		nym = mxGetScalar( mxGetField( prhs[0], 0, "nym" ) );
		}
	// nbx
	if(mxGetField( prhs[0], 0, "nbx" )!=NULL)
		{
		set_nbx = true;
		nbx = mxGetScalar( mxGetField( prhs[0], 0, "nbx" ) );
		}
	// nbu
	if(mxGetField( prhs[0], 0, "nbu" )!=NULL)
		{
		set_nbu = true;
		nbu = mxGetScalar( mxGetField( prhs[0], 0, "nbu" ) );
		}
	// ng
	if(mxGetField( prhs[0], 0, "ng" )!=NULL)
		{
		set_ng = true;
		ng = mxGetScalar( mxGetField( prhs[0], 0, "ng" ) );
		}
	// Vul
	if(mxGetField( prhs[0], 0, "Vul" )!=NULL)
		{
		set_Vul = true;
		Vul = mxGetPr( mxGetField( prhs[0], 0, "Vul" ) );
		}
	// Vxl
	if(mxGetField( prhs[0], 0, "Vxl" )!=NULL)
		{
		set_Vxl = true;
		Vxl = mxGetPr( mxGetField( prhs[0], 0, "Vxl" ) );
		}
	// Vxm
	if(mxGetField( prhs[0], 0, "Vxm" )!=NULL)
		{
		set_Vxm = true;
		Vxm = mxGetPr( mxGetField( prhs[0], 0, "Vxm" ) );
		}
	// Wl
	if(mxGetField( prhs[0], 0, "Wl" )!=NULL)
		{
		set_Wl = true;
		Wl = mxGetPr( mxGetField( prhs[0], 0, "Wl" ) );
		}
	// Wm
	if(mxGetField( prhs[0], 0, "Wm" )!=NULL)
		{
		set_Wm = true;
		Wm = mxGetPr( mxGetField( prhs[0], 0, "Wm" ) );
		}
	// yrl
	if(mxGetField( prhs[0], 0, "yrl" )!=NULL)
		{
		set_yrl = true;
		yrl = mxGetPr( mxGetField( prhs[0], 0, "yrl" ) );
		}
	// yrm
	if(mxGetField( prhs[0], 0, "yrm" )!=NULL)
		{
		set_yrm = true;
		yrm = mxGetPr( mxGetField( prhs[0], 0, "yrm" ) );
		}
	// x0
	if(mxGetField( prhs[0], 0, "x0" )!=NULL)
		{
		set_x0 = true;
		x0 = mxGetPr( mxGetField( prhs[0], 0, "x0" ) );
		}
	// Jbx
	if(mxGetField( prhs[0], 0, "Jbx" )!=NULL)
		{
		set_Jbx = true;
		Jbx = mxGetPr( mxGetField( prhs[0], 0, "Jbx" ) );
		}
	// lbx
	if(mxGetField( prhs[0], 0, "lbx" )!=NULL)
		{
		set_lbx = true;
		lbx = mxGetPr( mxGetField( prhs[0], 0, "lbx" ) );
		}
	// ubx
	if(mxGetField( prhs[0], 0, "ubx" )!=NULL)
		{
		set_ubx = true;
		ubx = mxGetPr( mxGetField( prhs[0], 0, "ubx" ) );
		}
	// Jbu
	if(mxGetField( prhs[0], 0, "Jbu" )!=NULL)
		{
		set_Jbu = true;
		Jbu = mxGetPr( mxGetField( prhs[0], 0, "Jbu" ) );
		}
	// lbu
	if(mxGetField( prhs[0], 0, "lbu" )!=NULL)
		{
		set_lbu = true;
		lbu = mxGetPr( mxGetField( prhs[0], 0, "lbu" ) );
		}
	// ubu
	if(mxGetField( prhs[0], 0, "ubu" )!=NULL)
		{
		set_ubu = true;
		ubu = mxGetPr( mxGetField( prhs[0], 0, "ubu" ) );
		}
	// C
	if(mxGetField( prhs[0], 0, "C" )!=NULL)
		{
		set_C = true;
		C = mxGetPr( mxGetField( prhs[0], 0, "C" ) );
		}
	// D
	if(mxGetField( prhs[0], 0, "D" )!=NULL)
		{
		set_D = true;
		D = mxGetPr( mxGetField( prhs[0], 0, "D" ) );
		}
	// lg
	if(mxGetField( prhs[0], 0, "lg" )!=NULL)
		{
		set_lg = true;
		lg = mxGetPr( mxGetField( prhs[0], 0, "lg" ) );
		}
	// ug
	if(mxGetField( prhs[0], 0, "ug" )!=NULL)
		{
		set_ug = true;
		ug = mxGetPr( mxGetField( prhs[0], 0, "ug" ) );
		}


	// opts_struct

	int N;							bool set_param_scheme_N = false;
	char * nlp_solver;
	char * qp_solver;
	int qp_solver_N_pcond;			bool set_qp_solver_N_pcond = false;
	char *sim_solver;
	int sim_solver_num_stages;		bool set_sim_solver_num_stages = false;
	int sim_solver_num_steps;		bool set_sim_solver_num_steps = false;

	// param_scheme_NN
	if(mxGetField( prhs[1], 0, "param_scheme_N" )!=NULL)
		{
		set_param_scheme_N = true;
		N = mxGetScalar( mxGetField( prhs[1], 0, "param_scheme_N" ) );
		}
	else
		{
		mexPrintf("\nerror: ocp_create: param_scheme_N not set!\n");
		return;
		}
	// nlp_solver
	// TODO check
	nlp_solver = mxArrayToString( mxGetField( prhs[1], 0, "nlp_solver" ) );
	// qp_solver
	// TODO check
	qp_solver = mxArrayToString( mxGetField( prhs[1], 0, "qp_solver" ) );
	// N_part_cond
	if(mxGetField( prhs[1], 0, "qp_solver_N_pcond" )!=NULL)
		{
		set_qp_solver_N_pcond = true;
		qp_solver_N_pcond = mxGetScalar( mxGetField( prhs[1], 0, "qp_solver_N_pcond" ) );
		}
	// sim_solver
	// TODO check
	sim_solver = mxArrayToString( mxGetField( prhs[1], 0, "sim_solver" ) );
	// sim_solver_num_stages
	if(mxGetField( prhs[1], 0, "sim_solver_num_stages" )!=NULL)
		{
		set_sim_solver_num_stages = true;
		sim_solver_num_stages = mxGetScalar( mxGetField( prhs[1], 0, "sim_solver_num_stages" ) );
		}
	// sim_solver_num_steps
	if(mxGetField( prhs[1], 0, "sim_solver_num_steps" )!=NULL)
		{
		set_sim_solver_num_steps = true;
		sim_solver_num_steps = mxGetScalar( mxGetField( prhs[1], 0, "sim_solver_num_steps" ) );
		}



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
	if(!strcmp(nlp_solver, "sqp"))
		{
		plan->nlp_solver = SQP;
		}
	else if(!strcmp(nlp_solver, "sqp_rti"))
		{
		plan->nlp_solver = SQP_RTI;
		}
	else
		{
		mexPrintf("\nnlp_solver not supported %s\n", nlp_solver);
		return;
		}
	
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
		if(!strcmp(sim_solver, "erk"))
			{
			plan->sim_solver_plan[ii].sim_solver = ERK;
			}
		else if(!strcmp(sim_solver, "irk"))
			{
			plan->sim_solver_plan[ii].sim_solver = IRK;
			}
		else
			{
			mexPrintf("\nsim_scheme not supported %s\n", sim_solver);
			return;
			}
//		plan->sim_solver_plan[ii].sim_solver = LIFTED_IRK;
//		plan->sim_solver_plan[ii].sim_solver = GNSF;
		}
	
	// constraints
	for(ii=0; ii<=N; ii++)
		{
		plan->nlp_constraints[ii] = BGH;
//		plan->nlp_constraints[ii] = BGHP;
		}

	// qp solver
	if(!strcmp(qp_solver, "partial_condensing_hpipm"))
		{
		plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
		}
	else if(!strcmp(qp_solver, "full_condensing_hpipm"))
		{
		plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_HPIPM;
		}
	else
		{
		mexPrintf("\nqp_solver not supported %s\n", qp_solver);
		return;
		}

	// TODO checks on initialization of plan !!!!!!!!

    ocp_nlp_config *config = ocp_nlp_config_create(*plan);


	ocp_nlp_plan_destroy(plan);
	


	/* dims */

	ocp_nlp_dims *dims = ocp_nlp_dims_create(config);
	// allocate tmp
	i_ptr = (int *) malloc((N+1)*sizeof(int));
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
	// nyl
	if(set_nyl)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dims_set_cost(config, dims, ii, "ny", &nyl);
			}
		}
	// nym
	if(set_nym)
		{
		ocp_nlp_dims_set_cost(config, dims, N, "ny", &nym);
		}
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
	// ng
	if(set_ng)
		{
		for(ii=0; ii<N; ii++)
			ocp_nlp_dims_set_constraints(config, dims, ii, "ng", &ng);
		}
			


	/* opts */

	void *opts = ocp_nlp_opts_create(config, dims);

	// qp_solver_N_pcond
	if(set_qp_solver_N_pcond)
		{
		ocp_nlp_opts_set(config, opts, "pcond_N2", &qp_solver_N_pcond);
		}
	// sim_solver_num_stages
	if(set_sim_solver_num_stages)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dynamics_opts_set(config, opts, ii, "num_stages", &sim_solver_num_stages);
			}
		}
	// sim_solver_num_steps
	if(set_sim_solver_num_steps)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_dynamics_opts_set(config, opts, ii, "num_steps", &sim_solver_num_steps);
			}
		}




	/* in */

	ocp_nlp_in *in = ocp_nlp_in_create(config, dims);

	// shooting nodes
	double Ts = T/N;
	for(ii=0; ii<N; ii++)
		{
		// TODO make setter for this
		in->Ts[ii] = Ts;
		}

	// cost: ls

	// lagrange term
	for(ii=0; ii<N; ii++)
		{
		if(set_Vul)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "Vu", Vul);
			}
		if(set_Vxl)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "Vx", Vxl);
			}
		if(set_Wl)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "W", Wl);
			}
		if(set_yrl)
			{
			ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yrl);
			}
		}
	// mayer term
	if(set_Vxm)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "Vx", Vxm);
		}
	if(set_Wm)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "W", Wm);
		}
	if(set_yrm)
		{
		ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yrm);
		}
	
	// constraints: b

	if(set_x0)
		{
		i_ptr = malloc(nx*sizeof(int));
		for(ii=0; ii<nx; ii++)
			{
			i_ptr[ii] = ii;
			}
		ocp_nlp_constraints_model_set(config, dims, in, 0, "idxbx", i_ptr);
		free(i_ptr);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", x0);
		ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", x0);
		}
	if(set_Jbx)
		{
		i_ptr = malloc(nbx*sizeof(int));
		for(ii=0; ii<nbx; ii++)
			{
			idx = -1;
			for(jj=0; jj<nx; jj++)
				{
				if(Jbx[ii+nbx*jj]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		for(ii=1; ii<=N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbx", i_ptr);
			}
		free(i_ptr);
		}
	if(set_lbx)
		{
		for(ii=1; ii<=N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "lbx", lbx);
			}
		}
	if(set_ubx)
		{
		for(ii=1; ii<=N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "ubx", ubx);
			}
		}
	if(set_Jbu)
		{
		i_ptr = malloc(nbu*sizeof(int));
		for(ii=0; ii<nbu; ii++)
			{
			idx = -1;
			for(jj=0; jj<nu; jj++)
				{
				if(Jbu[ii+nbu*jj]!=0.0)
					{
					i_ptr[ii] = jj;
					idx = jj;
					}
				}
			}
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbu", i_ptr);
			}
		free(i_ptr);
		}
	if(set_lbu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "lbu", lbu);
			}
		}
	if(set_ubu)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "ubu", ubu);
			}
		}
	if(set_C)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "C", C);
			}
		}
	if(set_D)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "D", D);
			}
		}
	if(set_lg)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "lg", lg);
			}
		}
	if(set_ug)
		{
		for(ii=0; ii<N; ii++)
			{
			ocp_nlp_constraints_model_set(config, dims, in, ii, "ug", ug);
			}
		}


#if 0
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
#endif


	/* out */

	ocp_nlp_out *out = ocp_nlp_out_create(config, dims);


	/* solver */

	ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, opts);



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


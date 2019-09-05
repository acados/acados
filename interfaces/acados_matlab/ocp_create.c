/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// acados
#include "acados_c/ocp_nlp_interface.h"
// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    long long *l_ptr;
    int *i_ptr;
    int *tmp_idx;
    double *d_ptr;
    char *c_ptr;

    int N, ii, jj, idx;
    char fun_name[50] = "ocp_create";
    char buffer [300]; // for error messages


    /* RHS */
    const mxArray *matlab_model = prhs[0];
    const mxArray *matlab_opts = prhs[1];


    /* LHS */
    // field names of output struct
    char *fieldnames[8];
    fieldnames[0] = (char*) mxMalloc(50);
    fieldnames[1] = (char*) mxMalloc(50);
    fieldnames[2] = (char*) mxMalloc(50);
    fieldnames[3] = (char*) mxMalloc(50);
    fieldnames[4] = (char*) mxMalloc(50);
    fieldnames[5] = (char*) mxMalloc(50);
    fieldnames[6] = (char*) mxMalloc(50);
    fieldnames[7] = (char*) mxMalloc(50);

    memcpy(fieldnames[0],"config",sizeof("config"));
    memcpy(fieldnames[1],"dims",sizeof("dims"));
    memcpy(fieldnames[2],"opts",sizeof("opts"));
    memcpy(fieldnames[3],"in",sizeof("in"));
    memcpy(fieldnames[4],"out",sizeof("out"));
    memcpy(fieldnames[5],"solver",sizeof("solver"));
    memcpy(fieldnames[6],"sens_out",sizeof("sens_out"));
    memcpy(fieldnames[7],"plan",sizeof("plan"));

    // create output struct
    plhs[0] = mxCreateStructMatrix(1, 1, 8, (const char **) fieldnames);

    mxFree( fieldnames[0] );
    mxFree( fieldnames[1] );
    mxFree( fieldnames[2] );
    mxFree( fieldnames[3] );
    mxFree( fieldnames[4] );
    mxFree( fieldnames[5] );
    mxFree( fieldnames[6] );
    mxFree( fieldnames[7] );


    /* plan & config */
    // param_scheme_N
    if (mxGetField( matlab_opts, 0, "param_scheme_N" )!=NULL)
    {
        N = mxGetScalar( mxGetField( matlab_opts, 0, "param_scheme_N" ) );
    }
    else
    {
       MEX_MISSING_ARGUMENT(fun_name, "param_scheme_N");
    }
    ocp_nlp_plan *plan = ocp_nlp_plan_create(N);


    // nlp solver
    char *nlp_solver = mxArrayToString( mxGetField( matlab_opts, 0, "nlp_solver" ) );
    if (!strcmp(nlp_solver, "sqp"))
    {
        plan->nlp_solver = SQP;
    }
    else if (!strcmp(nlp_solver, "sqp_rti"))
	{
        plan->nlp_solver = SQP_RTI;
	}
    else
	{
		MEX_SOLVER_NOT_SUPPORTED(fun_name, "nlp_solver", nlp_solver, "sqp, sqp_rti");
	}

    // cost type
	char *cost_type;
	ocp_nlp_cost_t cost_type_enum;
    if (mxGetField( matlab_model, 0, "cost_type" )!=NULL)
	{
		cost_type = mxArrayToString( mxGetField( matlab_model, 0, "cost_type" ) );
		if (!strcmp(cost_type, "linear_ls"))
		{
			cost_type_enum = LINEAR_LS;
		}
	    else if (!strcmp(cost_type, "nonlinear_ls"))
		{
			cost_type_enum = NONLINEAR_LS;
		}
	    else if (!strcmp(cost_type, "ext_cost"))
		{
			cost_type_enum = EXTERNALLY_PROVIDED;
		}
		else
		{
			MEX_SOLVER_NOT_SUPPORTED(fun_name, "cost_type", cost_type, "linear_ls, nonlinear_ls, ext_cost");
		}
	}
	else
	{
		MEX_MISSING_ARGUMENT(fun_name, "cost_type");
	}

	for (ii=0; ii<N; ii++)
	{
		plan->nlp_cost[ii] = cost_type_enum;
	}

    // cost e_type
	char *cost_type_e;
	ocp_nlp_cost_t cost_type_e_enum;
    if (mxGetField( matlab_model, 0, "cost_type_e" )!=NULL)
	{
		cost_type_e = mxArrayToString( mxGetField( matlab_model, 0, "cost_type_e" ) );
		if (!strcmp(cost_type_e, "linear_ls"))
		{
			cost_type_e_enum = LINEAR_LS;
		}
	    else if (!strcmp(cost_type_e, "nonlinear_ls"))
		{
			cost_type_e_enum = NONLINEAR_LS;
		}
	    else if (!strcmp(cost_type_e, "ext_cost"))
		{
			cost_type_e_enum = EXTERNALLY_PROVIDED;
		}
		else
		{
			MEX_SOLVER_NOT_SUPPORTED(fun_name, "cost_type_e", cost_type_e, "linear_ls, nonlinear_ls, ext_cost");
		}
	}
	else
	{
		MEX_MISSING_ARGUMENT(fun_name, "cost_type_e");
	}
	plan->nlp_cost[N] = cost_type_e_enum;



    // dynamics: dyn_type, sim_method
	char *dyn_type;
	ocp_nlp_dynamics_t dyn_type_enum;
    if (mxGetField( matlab_model, 0, "dyn_type" )!=NULL)
	{
        dyn_type = mxArrayToString( mxGetField( matlab_model, 0, "dyn_type" ) );
		if (!strcmp(dyn_type, "explicit"))
		{
			dyn_type_enum = CONTINUOUS_MODEL;
		}
		else if (!strcmp(dyn_type, "implicit"))
		{
			dyn_type_enum = CONTINUOUS_MODEL;
		}
	    else if (!strcmp(dyn_type, "discrete"))
		{
			dyn_type_enum = DISCRETE_MODEL;
		}
		else
		{
			MEX_SOLVER_NOT_SUPPORTED(fun_name, "dyn_type", dyn_type, "explicit, implicit, discrete");
		}
	}
	else
	{
		MEX_MISSING_ARGUMENT(fun_name, "dyn_type");
	}
	for (ii=0; ii<N; ii++)
	{
		plan->nlp_dynamics[ii] = dyn_type_enum;
	}

    char *sim_method;
	sim_solver_t sim_method_enum;

	if (strcmp(dyn_type, "discrete")!=0) // discrete dont need integrator
	{
		if (mxGetField( matlab_opts, 0, "sim_method" )!=NULL)
		{
			sim_method = mxArrayToString( mxGetField( matlab_opts, 0, "sim_method" ) );
		}
		else
		{
			MEX_MISSING_ARGUMENT(fun_name, "sim_method");
		}
	}
	else
	{
		sim_method_enum = INVALID_SIM_SOLVER; // discrete
	}

	if (!strcmp(dyn_type, "explicit"))
	{
		if (!strcmp(sim_method, "erk"))
		{
			sim_method_enum = ERK;
		}
		else
		{
			MEX_SOLVER_NOT_SUPPORTED_GIVEN(fun_name, "sim_method", sim_method, "explicit dynamics", "erk");
		}
	}
    else if (!strcmp(dyn_type, "implicit"))
	{
        if (!strcmp(sim_method, "irk"))
		{
			sim_method_enum = IRK;
		}
        else if (!strcmp(sim_method, "irk_gnsf"))
		{
			sim_method_enum = GNSF;
		}
        else
		{
			MEX_SOLVER_NOT_SUPPORTED_GIVEN(fun_name, "sim_method", sim_method, "implicit dynamics", "irk, irk_gnsf");
		}
	}
	for (ii=0; ii<N; ii++)
	{
	    plan->sim_solver_plan[ii].sim_solver = sim_method_enum;
	}


    // constraints
    char *constr_type;
    if (mxGetField( matlab_model, 0, "constr_type" )!=NULL)
	{
        constr_type = mxArrayToString( mxGetField( matlab_model, 0, "constr_type" ) );
	}
	else
	{
		MEX_MISSING_ARGUMENT(fun_name, "constr_type");
	}

    if (!strcmp(constr_type, "bgh"))
	{
        for (ii=0; ii<=N; ii++)
		{
            plan->nlp_constraints[ii] = BGH;
		}
	}
    else
	{
		MEX_SOLVER_NOT_SUPPORTED(fun_name, "constr_type", constr_type, "bgh");
	}

    // qp solver
	char *qp_solver;
    if (mxGetField( matlab_opts, 0, "qp_solver" )!=NULL)
	{
        qp_solver = mxArrayToString( mxGetField( matlab_opts, 0, "qp_solver" ) );
	}
	else
	{
		MEX_MISSING_ARGUMENT(fun_name, "qp_solver");
	}

    if (!strcmp(qp_solver, "partial_condensing_hpipm"))
	{
        plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
	}
    else if (!strcmp(qp_solver, "full_condensing_hpipm"))
	{
        plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_HPIPM;
	}
#if defined(ACADOS_WITH_QPOASES)
    else if (!strcmp(qp_solver, "full_condensing_qpoases"))
	{
        plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_QPOASES;
	}
#endif
    else
	{
		MEX_SOLVER_NOT_SUPPORTED(fun_name, "qp_solver", qp_solver,
		     "partial_condensing_hpipm, full_condensing_hpipm, full_condensing_qpoases");
	}
    
    // regularization
	char *regularize_method;
    if (mxGetField( matlab_opts, 0, "regularize_method" )!=NULL)
	{
        regularize_method = mxArrayToString( mxGetField( matlab_opts, 0, "regularize_method" ) );
	}
	else
	{
		MEX_MISSING_ARGUMENT(fun_name, "regularize_method");
	}

    if (!strcmp(regularize_method, "no_regularize"))
	{
        plan->regularization = NO_REGULARIZE;
	}
    else if (!strcmp(regularize_method, "mirror"))
	{
        plan->regularization = MIRROR;
	}
    else if (!strcmp(regularize_method, "project"))
	{
        plan->regularization = PROJECT;
	}
    else if (!strcmp(regularize_method, "project_reduc_hess"))
	{
        plan->regularization = PROJECT_REDUC_HESS;
	}
    else if (!strcmp(regularize_method, "convexify"))
	{
		plan->regularization = CONVEXIFY;
	}
    else
	{
		MEX_SOLVER_NOT_SUPPORTED(fun_name, "regularize_method", regularize_method,
		     "no_regularize, mirror, project, project_reduc_hess, convexify");
	}

    ocp_nlp_config *config = ocp_nlp_config_create(*plan);


    /* dims */
    int nx, nu;
    int nz = 0;
    int ny, ny_e;
    int nbx;
    int nbu;
    int ng;
    int ng_e;
    int nh;
    int nh_e;
    int ns = 0;
    int ns_e = 0;
    int nsbu = 0;
    int nsbx = 0;
    int nsg = 0;
    int nsg_e = 0;
    int nsh = 0;
    int nsh_e = 0;


    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);

    // nx
    if (mxGetField( matlab_model, 0, "dim_nx" )!=NULL)
	{
        nx = mxGetScalar( mxGetField( prhs[0], 0, "dim_nx" ) );
	}
    else
	{
        MEX_MISSING_ARGUMENT(fun_name, "dim_nx");
	}

    i_ptr = (int *) malloc((N+1)*sizeof(int));
    for (ii=0; ii<=N; ii++)
        i_ptr[ii] = nx;
    ocp_nlp_dims_set_opt_vars(config, dims, "nx", i_ptr);

    // nu
    if (mxGetField( matlab_model, 0, "dim_nu" )!=NULL)
	{
        nu = mxGetScalar( mxGetField( prhs[0], 0, "dim_nu" ) );
	}
    else
	{
        MEX_MISSING_ARGUMENT(fun_name, "dim_nu");
	}
    for (ii=0; ii<N; ii++)
        i_ptr[ii] = nu;
    i_ptr[N] = 0;
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", i_ptr);
	
	// nz
    if (mxGetField( prhs[0], 0, "dim_nz" )!=NULL)
	{
        nz = mxGetScalar( mxGetField( prhs[0], 0, "dim_nz" ) );
        for (ii=0; ii<=N; ii++)
            i_ptr[ii] = nz;
        ocp_nlp_dims_set_opt_vars(config, dims, "nz", i_ptr);
	}
    free(i_ptr);

    // gnsf stuff
    int gnsf_nx1, gnsf_nz1, gnsf_nuhat, gnsf_ny, gnsf_nout;
    if (!strcmp(sim_method, "irk_gnsf"))
	{
		// nx1
        if (mxGetField( matlab_model, 0, "dim_gnsf_nx1" )!=NULL)
		{
            gnsf_nx1 = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nx1" ) );
		}
		else
		{
			MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nx1", "irk_gnsf")
		}
		// nz1
        if (mxGetField( matlab_model, 0, "dim_gnsf_nz1" )!=NULL)
		{
            gnsf_nz1 = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nz1" ) );
		}
		else
		{
			MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nz1", "irk_gnsf")
		}
		// nuhat
		if (mxGetField( matlab_model, 0, "dim_gnsf_nuhat" )!=NULL)
		{
            gnsf_nuhat = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nuhat" ) );
		}
		else
		{
			MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nuhat", "irk_gnsf")
		}
		// ny
        if (mxGetField( matlab_model, 0, "dim_gnsf_ny" )!=NULL)
		{
            gnsf_ny = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_ny" ) );
		}
		else
		{
			MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_ny", "irk_gnsf")
		}
		// nout
        if (mxGetField( matlab_model, 0, "dim_gnsf_nout" )!=NULL)
		{
            gnsf_nout = mxGetScalar( mxGetField( matlab_model, 0, "dim_gnsf_nout" ) );
		}
		else
		{
			MEX_MISSING_ARGUMENT_MODULE(fun_name, "dim_gnsf_nout", "irk_gnsf")
		}
		// assign
		for (ii=0; ii<N; ii++)
		{
			ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nx1", &gnsf_nx1);
			ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nz1", &gnsf_nz1);
			ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nuhat", &gnsf_nuhat);
			ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_ny", &gnsf_ny);
			ocp_nlp_dims_set_dynamics(config, dims, ii, "gnsf_nout", &gnsf_nout);
		}
	}


    // ny
    if (mxGetField( matlab_model, 0, "dim_ny" )!=NULL)
	{
        ny = mxGetScalar( mxGetField( matlab_model, 0, "dim_ny" ) );
		for (ii=0; ii<N; ii++)
		{
            ocp_nlp_dims_set_cost(config, dims, ii, "ny", &ny);
		}
	}
    // ny_e
    if (mxGetField( matlab_model, 0, "dim_ny_e" )!=NULL)
	{
        ny_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_ny_e" ) );
        ocp_nlp_dims_set_cost(config, dims, N, "ny", &ny_e);
	}

	// constraint dims
	// nbx
    ocp_nlp_dims_set_constraints(config, dims, 0, "nbx", &nx);
    if (mxGetField( matlab_model, 0, "dim_nbx" )!=NULL)
	{
        nbx = mxGetScalar( mxGetField( matlab_model, 0, "dim_nbx" ) );
        for (ii=1; ii<=N; ii++)
		{
            ocp_nlp_dims_set_constraints(config, dims, ii, "nbx", &nbx);
		}
	}
    // nbu
    if (mxGetField( matlab_model, 0, "dim_nbu" )!=NULL)
	{
        nbu = mxGetScalar( mxGetField( matlab_model, 0, "dim_nbu" ) );
        for (ii=0; ii<N; ii++)
		{
            ocp_nlp_dims_set_constraints(config, dims, ii, "nbu", &nbu);
		}
	}
    // ng
    if (mxGetField( matlab_model, 0, "dim_ng" )!=NULL)
	{
        ng = mxGetScalar( mxGetField( matlab_model, 0, "dim_ng" ) );
        for (ii=0; ii<N; ii++)
		{
            ocp_nlp_dims_set_constraints(config, dims, ii, "ng", &ng);
        }
	}
    // ng_e
    if (mxGetField( matlab_model, 0, "dim_ng_e" )!=NULL)
	{
        ng_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_ng_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "ng", &ng_e);
	}

    // nh
    if (mxGetField( matlab_model, 0, "dim_nh" )!=NULL)
	{
        nh = mxGetScalar( mxGetField( matlab_model, 0, "dim_nh" ) );
        for (ii=0; ii<N; ii++)
		{
            ocp_nlp_dims_set_constraints(config, dims, ii, "nh", &nh);
        }
	}
    // nh_e
    if (mxGetField( matlab_model, 0, "dim_nh_e" )!=NULL)
	{
        nh_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nh_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nh", &nh_e);
	}

	// slack dims
    // ns
    if (mxGetField( matlab_model, 0, "dim_ns" )!=NULL)
	{
        ns = mxGetScalar( mxGetField( matlab_model, 0, "dim_ns" ) );
	}

    // ns_e
    if (mxGetField( matlab_model, 0, "dim_ns_e" )!=NULL)
	{
        ns_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_ns_e" ) );
	}
    // nsbu
    if (mxGetField( matlab_model, 0, "dim_nsbu" )!=NULL)
	{
        nsbu = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsbu" ) );
		for (ii=0; ii<N; ii++)
		{
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsbu", &nsbu);
		}
	}
    // nsbx
    if (mxGetField( matlab_model, 0, "dim_nsbx" )!=NULL)
	{
        nsbx = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsbx" ) );
        for (ii=1; ii<=N; ii++) // TODO fix stage 0 !!!!!
		{
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsbx", &nsbx);
		}
	}
    // nsg
    if (mxGetField( matlab_model, 0, "dim_nsg" )!=NULL)
	{
        nsg = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsg" ) );
        for (ii=0; ii<N; ii++)
		{
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsg", &nsg);
		}
	}
    // nsg_e
    if (mxGetField( matlab_model, 0, "dim_nsg_e" )!=NULL)
	{
        nsg_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsg_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nsg", &nsg_e);
	}
    // nsh
    if (mxGetField( matlab_model, 0, "dim_nsh" )!=NULL)
	{
        nsh = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsh" ) );
        for (ii=0; ii<N; ii++)
        {
            ocp_nlp_dims_set_constraints(config, dims, ii, "nsh", &nsh);
        }
	}
    // nsh_e
    if (mxGetField( matlab_model, 0, "dim_nsh_e" )!=NULL)
	{
        nsh_e = mxGetScalar( mxGetField( matlab_model, 0, "dim_nsh_e" ) );
        ocp_nlp_dims_set_constraints(config, dims, N, "nsh", &nsh_e);
	}
    // ns
    if (ns!=nsbu+nsbx+nsg+nsh)
	{
		sprintf(buffer,"ocp_create: ns!=nsbu+nsbx+nsg+nsh, got ns=%d, nsbu=%d, nsbx=%d nsg=%d nsh=%d\n",
			ns, nsbu, nsbx, nsg, nsh);
        mexErrMsgTxt(buffer);
	}
    if (ns_e!=nsbx+nsg_e+nsh_e)
	{
		sprintf(buffer,"ocp_create: ns_e!=nsbx+nsg_e+nsh_e, got ns_e=%d, nsbx=%d nsg_e=%d nsh_e=%d\n",
			ns, nsbx, nsg_e, nsh_e);
        mexErrMsgTxt(buffer);
	}

	i_ptr = (int *) malloc((N+1)*sizeof(int));
    // TODO fix stage 0 !!!!!!!!
    i_ptr[0] = nsbu+nsg+nsh; // XXX not nsbx !!!!!!!!!!
    for (ii=1; ii<N; ii++)
        i_ptr[ii] = ns;
    i_ptr[N] = ns_e;
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", i_ptr);
	free(i_ptr);

    /* opts */

    void *opts = ocp_nlp_opts_create(config, dims);

    int nlp_solver_max_iter;
    int nlp_solver_ext_qp_res;                bool set_nlp_solver_ext_qp_res = false;
    int qp_solver_iter_max;                    bool set_qp_solver_iter_max = false;
    int qp_solver_cond_N;                    bool set_qp_solver_cond_N = false;
    int qp_solver_cond_ric_alg;                bool set_qp_solver_cond_ric_alg = false;
    int qp_solver_ric_alg;                    bool set_qp_solver_ric_alg = false;
    int qp_solver_warm_start;                bool set_qp_solver_warm_start = false;
    int sim_method_num_stages;                bool set_sim_method_num_stages = false;
    int sim_method_num_steps;                bool set_sim_method_num_steps = false;


    // nlp solver exact hessian
    bool nlp_solver_exact_hessian = false;
    c_ptr = mxArrayToString( mxGetField( matlab_opts, 0, "nlp_solver_exact_hessian" ) );
	MEX_STR_TO_BOOL(fun_name, c_ptr, &nlp_solver_exact_hessian, "nlp_solver_exact_hessian");
	// TODO: check if possible: e.g. not with irk_gnsf
	ocp_nlp_opts_set(config, opts, "exact_hess", &nlp_solver_exact_hessian);


    // nlp solver max iter
    if (mxGetField( matlab_opts, 0, "nlp_solver_max_iter" )!=NULL)
	{
        nlp_solver_max_iter = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_max_iter" ) );
        ocp_nlp_opts_set(config, opts, "max_iter", &nlp_solver_max_iter);
	}
    // nlp solver exit tolerances
    double nlp_solver_tol_stat, nlp_solver_tol_eq, nlp_solver_tol_ineq, nlp_solver_tol_comp;
    if (mxGetField( matlab_opts, 0, "nlp_solver_tol_stat" )!=NULL)
	{
        nlp_solver_tol_stat = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_stat" ) );
        ocp_nlp_opts_set(config, opts, "tol_stat", &nlp_solver_tol_stat);
	}
    if (mxGetField( matlab_opts, 0, "nlp_solver_tol_eq" )!=NULL)
	{
        nlp_solver_tol_eq = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_eq" ) );
        ocp_nlp_opts_set(config, opts, "tol_eq", &nlp_solver_tol_eq);
	}
    if (mxGetField( matlab_opts, 0, "nlp_solver_tol_ineq" )!=NULL)
	{
        nlp_solver_tol_ineq = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_ineq" ) );
        ocp_nlp_opts_set(config, opts, "tol_ineq", &nlp_solver_tol_ineq);
	}
    if (mxGetField( matlab_opts, 0, "nlp_solver_tol_comp" )!=NULL)
	{
        nlp_solver_tol_comp = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_tol_comp" ) );
        ocp_nlp_opts_set(config, opts, "tol_comp", &nlp_solver_tol_comp);
	}


    // nlp solver ext qp res
    if (mxGetField( matlab_opts, 0, "nlp_solver_ext_qp_res" )!=NULL)
        {
        set_nlp_solver_ext_qp_res = true;
        nlp_solver_ext_qp_res = mxGetScalar( mxGetField( matlab_opts, 0, "nlp_solver_ext_qp_res" ) );
        }
    // qp_solver
    qp_solver = mxArrayToString( mxGetField( matlab_opts, 0, "qp_solver" ) );
    // iter_max
    if (mxGetField( matlab_opts, 0, "qp_solver_iter_max" )!=NULL)
        {
        set_qp_solver_iter_max = true;
        qp_solver_iter_max = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_iter_max" ) );
        }
    // N_part_cond
    if (mxGetField( matlab_opts, 0, "qp_solver_cond_N" )!=NULL)
        {
        set_qp_solver_cond_N = true;
        qp_solver_cond_N = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_cond_N" ) );
        }
    // cond riccati-like algorithm
    if (mxGetField( matlab_opts, 0, "qp_solver_cond_ric_alg" )!=NULL)
        {
        set_qp_solver_cond_ric_alg = true;
        qp_solver_cond_ric_alg = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_cond_ric_alg" ) );
        }
    // hpipm: riccati algorithm
    if (mxGetField( matlab_opts, 0, "qp_solver_ric_alg" )!=NULL)
        {
        set_qp_solver_ric_alg = true;
        qp_solver_ric_alg = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_ric_alg" ) );
        }
    // qp solver: warm start
    if (mxGetField( matlab_opts, 0, "qp_solver_warm_start" )!=NULL)
        {
        set_qp_solver_warm_start = true;
        qp_solver_warm_start = mxGetScalar( mxGetField( matlab_opts, 0, "qp_solver_warm_start" ) );
        }
    // sim_method_num_stages
    if (mxGetField( matlab_opts, 0, "sim_method_num_stages" )!=NULL)
        {
        set_sim_method_num_stages = true;
        sim_method_num_stages = mxGetScalar( mxGetField( matlab_opts, 0, "sim_method_num_stages" ) );
        }
    // sim_method_num_steps
    if (mxGetField( matlab_opts, 0, "sim_method_num_steps" )!=NULL)
        {
        set_sim_method_num_steps = true;
        sim_method_num_steps = mxGetScalar( mxGetField( matlab_opts, 0, "sim_method_num_steps" ) );
        }




    // model
    double T;        bool set_T = false;
    // cost
    double *Vu;        bool set_Vu = false;
    double *Vx;        bool set_Vx = false;
    double *Vx_e;    bool set_Vx_e = false;
    double *W;        bool set_W = false;
    double *W_e;    bool set_W_e = false;
    double *yr;        bool set_yr = false;
    double *yr_e;    bool set_yr_e = false;
    double *Z;        bool set_Z = false;
    double *Z_e;    bool set_Z_e = false;
    double *Zl;        bool set_Zl = false;
    double *Zl_e;    bool set_Zl_e = false;
    double *Zu;        bool set_Zu = false;
    double *Zu_e;    bool set_Zu_e = false;
    double *z;        bool set_z = false;
    double *z_e;    bool set_z_e = false;
    double *zl;        bool set_zl = false;
    double *zl_e;    bool set_zl_e = false;
    double *zu;        bool set_zu = false;
    double *zu_e;    bool set_zu_e = false;
    // constraints
    double *x0;        bool set_x0 = false;
    double *Jbx;    bool set_Jbx = false;
    double *lbx;    bool set_lbx = false;
    double *ubx;    bool set_ubx = false;
    double *Jbu;    bool set_Jbu = false;
    double *lbu;    bool set_lbu = false;
    double *ubu;    bool set_ubu = false;
    double *C;        bool set_C = false;
    double *D;        bool set_D = false;
    double *lg;        bool set_lg = false;
    double *ug;        bool set_ug = false;
    double *C_e;    bool set_C_e = false;
    double *lg_e;    bool set_lg_e = false;
    double *ug_e;    bool set_ug_e = false;
    double *lh;        bool set_lh = false;
    double *uh;        bool set_uh = false;
    double *lh_e;    bool set_lh_e = false;
    double *uh_e;    bool set_uh_e = false;
    double *Jsbu;    bool set_Jsbu = false;
//    double *lsbu;    bool set_lsbu = false;
//    double *usbu;    bool set_usbu = false;
    double *Jsbx;    bool set_Jsbx = false;
//    double *lsbx;    bool set_lsbx = false;
//    double *usbx;    bool set_usbx = false;
    double *Jsg;    bool set_Jsg = false;
//    double *lsg;    bool set_lsg = false;
//    double *usg;    bool set_usg = false;
    double *Jsg_e;    bool set_Jsg_e = false;
//    double *lsg_e;    bool set_lsg_e = false;
//    double *usg_e;    bool set_usg_e = false;
    double *Jsh;    bool set_Jsh = false;
//    double *lsh;    bool set_lsh = false;
//    double *ush;    bool set_ush = false;
    double *Jsh_e;    bool set_Jsh_e = false;
//    double *lsh_e;    bool set_lsh_e = false;
//    double *ush_e;    bool set_ush_e = false;
    // dynamics
    // trajectory initialization
    double *x_init; bool set_x_init = false;
    double *u_init; bool set_u_init = false;

    // dims
    // T
    if (mxGetField( matlab_model, 0, "T" )!=NULL)
        {
        set_T = true;
        T = mxGetScalar( mxGetField( matlab_model, 0, "T" ) );
        }
    else
        {
        mexPrintf("\nerror: ocp_create: T not set!\n");
        return;
        }










    // Vu
    if (mxGetField( matlab_model, 0, "cost_Vu" )!=NULL)
        {
        set_Vu = true;
        Vu = mxGetPr( mxGetField( matlab_model, 0, "cost_Vu" ) );
        }
    // Vx
    if (mxGetField( matlab_model, 0, "cost_Vx" )!=NULL)
        {
        set_Vx = true;
        Vx = mxGetPr( mxGetField( matlab_model, 0, "cost_Vx" ) );
        }
    // Vx_e
    if (mxGetField( matlab_model, 0, "cost_Vx_e" )!=NULL)
        {
        set_Vx_e = true;
        Vx_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Vx_e" ) );
        }
    // W
    if (mxGetField( matlab_model, 0, "cost_W" )!=NULL)
        {
        set_W = true;
        W = mxGetPr( mxGetField( matlab_model, 0, "cost_W" ) );
        }
    // W_e
    if (mxGetField( matlab_model, 0, "cost_W_e" )!=NULL)
        {
        set_W_e = true;
        W_e = mxGetPr( mxGetField( matlab_model, 0, "cost_W_e" ) );
        }
    // yr
    if (mxGetField( matlab_model, 0, "cost_yr" )!=NULL)
        {
        set_yr = true;
        yr = mxGetPr( mxGetField( matlab_model, 0, "cost_yr" ) );
        }
    // yr_e
    if (mxGetField( matlab_model, 0, "cost_yr_e" )!=NULL)
        {
        set_yr_e = true;
        yr_e = mxGetPr( mxGetField( matlab_model, 0, "cost_yr_e" ) );
        }
    // Z
    if (mxGetField( matlab_model, 0, "cost_Z" )!=NULL)
        {
        set_Z = true;
        Z = mxGetPr( mxGetField( matlab_model, 0, "cost_Z" ) );
        }
    // Z_e
    if (mxGetField( matlab_model, 0, "cost_Z_e" )!=NULL)
        {
        set_Z_e = true;
        Z_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Z_e" ) );
        }
    // Zl
    if (mxGetField( matlab_model, 0, "cost_Zl" )!=NULL)
        {
        set_Zl = true;
        Zl = mxGetPr( mxGetField( matlab_model, 0, "cost_Zl" ) );
        }
    // Zl_e
    if (mxGetField( matlab_model, 0, "cost_Zl_e" )!=NULL)
        {
        set_Zl_e = true;
        Zl_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Zl_e" ) );
        }
    // Zu
    if (mxGetField( matlab_model, 0, "cost_Zu" )!=NULL)
        {
        set_Zu = true;
        Zu = mxGetPr( mxGetField( matlab_model, 0, "cost_Zu" ) );
        }
    // Zu_e
    if (mxGetField( matlab_model, 0, "cost_Zu_e" )!=NULL)
        {
        set_Zu_e = true;
        Zu_e = mxGetPr( mxGetField( matlab_model, 0, "cost_Zu_e" ) );
        }
    // z
    if (mxGetField( matlab_model, 0, "cost_z" )!=NULL)
        {
        set_z = true;
        z = mxGetPr( mxGetField( matlab_model, 0, "cost_z" ) );
        }
    // z_e
    if (mxGetField( matlab_model, 0, "cost_z_e" )!=NULL)
        {
        set_z_e = true;
        z_e = mxGetPr( mxGetField( matlab_model, 0, "cost_z_e" ) );
        }
    // zl
    if (mxGetField( matlab_model, 0, "cost_zl" )!=NULL)
        {
        set_zl = true;
        zl = mxGetPr( mxGetField( matlab_model, 0, "cost_zl" ) );
        }
    // zl_e
    if (mxGetField( matlab_model, 0, "cost_zl_e" )!=NULL)
        {
        set_zl_e = true;
        zl_e = mxGetPr( mxGetField( matlab_model, 0, "cost_zl_e" ) );
        }
    // zu
    if (mxGetField( matlab_model, 0, "cost_zu" )!=NULL)
        {
        set_zu = true;
        zu = mxGetPr( mxGetField( matlab_model, 0, "cost_zu" ) );
        }
    // zu_e
    if (mxGetField( matlab_model, 0, "cost_zu_e" )!=NULL)
        {
        set_zu_e = true;
        zu_e = mxGetPr( mxGetField( matlab_model, 0, "cost_zu_e" ) );
        }
    // x0
    if (mxGetField( matlab_model, 0, "constr_x0" )!=NULL)
        {
        set_x0 = true;
        x0 = mxGetPr( mxGetField( matlab_model, 0, "constr_x0" ) );
        }
    // Jbx
    if (mxGetField( matlab_model, 0, "constr_Jbx" )!=NULL)
        {
        set_Jbx = true;
        Jbx = mxGetPr( mxGetField( matlab_model, 0, "constr_Jbx" ) );
        }
    // lbx
    if (mxGetField( matlab_model, 0, "constr_lbx" )!=NULL)
        {
        set_lbx = true;
        lbx = mxGetPr( mxGetField( matlab_model, 0, "constr_lbx" ) );
        }
    // ubx
    if (mxGetField( matlab_model, 0, "constr_ubx" )!=NULL)
        {
        set_ubx = true;
        ubx = mxGetPr( mxGetField( matlab_model, 0, "constr_ubx" ) );
        }
    // Jbu
    if (mxGetField( matlab_model, 0, "constr_Jbu" )!=NULL)
        {
        set_Jbu = true;
        Jbu = mxGetPr( mxGetField( matlab_model, 0, "constr_Jbu" ) );
        }
    // lbu
    if (mxGetField( matlab_model, 0, "constr_lbu" )!=NULL)
        {
        set_lbu = true;
        lbu = mxGetPr( mxGetField( matlab_model, 0, "constr_lbu" ) );
        }
    // ubu
    if (mxGetField( matlab_model, 0, "constr_ubu" )!=NULL)
        {
        set_ubu = true;
        ubu = mxGetPr( mxGetField( matlab_model, 0, "constr_ubu" ) );
        }
    // C
    if (mxGetField( matlab_model, 0, "constr_C" )!=NULL)
        {
        set_C = true;
        C = mxGetPr( mxGetField( matlab_model, 0, "constr_C" ) );
        }
    // D
    if (mxGetField( matlab_model, 0, "constr_D" )!=NULL)
        {
        set_D = true;
        D = mxGetPr( mxGetField( matlab_model, 0, "constr_D" ) );
        }
    // lg
    if (mxGetField( matlab_model, 0, "constr_lg" )!=NULL)
        {
        set_lg = true;
        lg = mxGetPr( mxGetField( matlab_model, 0, "constr_lg" ) );
        }
    // ug
    if (mxGetField( matlab_model, 0, "constr_ug" )!=NULL)
        {
        set_ug = true;
        ug = mxGetPr( mxGetField( matlab_model, 0, "constr_ug" ) );
        }
    // C_e
    if (mxGetField( matlab_model, 0, "constr_C_e" )!=NULL)
        {
        set_C_e = true;
        C_e = mxGetPr( mxGetField( matlab_model, 0, "constr_C_e" ) );
        }
    // lg_e
    if (mxGetField( matlab_model, 0, "constr_lg_e" )!=NULL)
        {
        set_lg_e = true;
        lg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lg_e" ) );
        }
    // ug_e
    if (mxGetField( matlab_model, 0, "constr_ug_e" )!=NULL)
        {
        set_ug_e = true;
        ug_e = mxGetPr( mxGetField( matlab_model, 0, "constr_ug_e" ) );
        }
    // lh
    if (mxGetField( matlab_model, 0, "constr_lh" )!=NULL)
        {
        set_lh = true;
        lh = mxGetPr( mxGetField( matlab_model, 0, "constr_lh" ) );
        }
    // uh
    if (mxGetField( matlab_model, 0, "constr_uh" )!=NULL)
        {
        set_uh = true;
        uh = mxGetPr( mxGetField( matlab_model, 0, "constr_uh" ) );
        }
    // lh_e
    if (mxGetField( matlab_model, 0, "constr_lh_e" )!=NULL)
        {
        set_lh_e = true;
        lh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lh_e" ) );
        }
    // uh_e
    if (mxGetField( matlab_model, 0, "constr_uh_e" )!=NULL)
        {
        set_uh_e = true;
        uh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_uh_e" ) );
        }
    // Jsbu
    if (mxGetField( matlab_model, 0, "constr_Jsbu" )!=NULL)
        {
        set_Jsbu = true;
        Jsbu = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsbu" ) );
        }
    // lsbu
//    if (mxGetField( matlab_model, 0, "constr_lsbu" )!=NULL)
//        {
//        set_lsbu = true;
//        lsbu = mxGetPr( mxGetField( matlab_model, 0, "constr_lsbu" ) );
//        }
    // usbu
//    if (mxGetField( matlab_model, 0, "constr_usbu" )!=NULL)
//        {
//        set_usbu = true;
//        usbu = mxGetPr( mxGetField( matlab_model, 0, "constr_usbu" ) );
//        }
    // Jsbx
    if (mxGetField( matlab_model, 0, "constr_Jsbx" )!=NULL)
        {
        set_Jsbx = true;
        Jsbx = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsbx" ) );
        }
    // lsbx
//    if (mxGetField( matlab_model, 0, "constr_lsbx" )!=NULL)
//        {
//        set_lsbx = true;
//        lsbx = mxGetPr( mxGetField( matlab_model, 0, "constr_lsbx" ) );
//        }
    // usbx
//    if (mxGetField( matlab_model, 0, "constr_usbx" )!=NULL)
//        {
//        set_usbx = true;
//        usbx = mxGetPr( mxGetField( matlab_model, 0, "constr_usbx" ) );
//        }
    // Jsg
    if (mxGetField( matlab_model, 0, "constr_Jsg" )!=NULL)
        {
        set_Jsg = true;
        Jsg = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsg" ) );
        }
    // lsg
//    if (mxGetField( matlab_model, 0, "constr_lsg" )!=NULL)
//        {
//        set_lsg = true;
//        lsg = mxGetPr( mxGetField( matlab_model, 0, "constr_lsg" ) );
//        }
    // usg
//    if (mxGetField( matlab_model, 0, "constr_usg" )!=NULL)
//        {
//        set_usg = true;
//        usg = mxGetPr( mxGetField( matlab_model, 0, "constr_usg" ) );
//        }
    // Jsg_e
    if (mxGetField( matlab_model, 0, "constr_Jsg_e" )!=NULL)
        {
        set_Jsg_e = true;
        Jsg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsg_e" ) );
        }
    // lsg_e
//    if (mxGetField( matlab_model, 0, "constr_lsg_e" )!=NULL)
//        {
//        set_lsg_e = true;
//        lsg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lsg_e" ) );
//        }
    // usg_e
//    if (mxGetField( matlab_model, 0, "constr_usg_e" )!=NULL)
//        {
//        set_usg_e = true;
//        usg_e = mxGetPr( mxGetField( matlab_model, 0, "constr_usg_e" ) );
//        }
    // Jsh
    if (mxGetField( matlab_model, 0, "constr_Jsh" )!=NULL)
        {
        set_Jsh = true;
        Jsh = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsh" ) );
        }
    // lsh
//    if (mxGetField( matlab_model, 0, "constr_lsh" )!=NULL)
//        {
//        set_lsh = true;
//        lsh = mxGetPr( mxGetField( matlab_model, 0, "constr_lsh" ) );
//        }
    // ush
//    if (mxGetField( matlab_model, 0, "constr_ush" )!=NULL)
//        {
//        set_ush = true;
//        ush = mxGetPr( mxGetField( matlab_model, 0, "constr_ush" ) );
//        }
    // Jsh_e
    if (mxGetField( matlab_model, 0, "constr_Jsh_e" )!=NULL)
        {
        set_Jsh_e = true;
        Jsh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_Jsh_e" ) );
        }
    // lsh_e
//    if (mxGetField( matlab_model, 0, "constr_lsh_e" )!=NULL)
//        {
//        set_lsh_e = true;
//        lsh_e = mxGetPr( mxGetField( matlab_model, 0, "constr_lsh_e" ) );
//        }
    // ush_e
//    if (mxGetField( matlab_model, 0, "constr_ush_e" )!=NULL)
//        {
//        set_ush_e = true;
//        ush_e = mxGetPr( mxGetField( matlab_model, 0, "constr_ush_e" ) );
//        }

    // trajectory initialization
    // u_init
    if (mxGetField( matlab_model, 0, "u_init" )!=NULL)
        {
        set_u_init = true;
        u_init = mxGetPr( mxGetField( matlab_model, 0, "u_init" ) );
        }
    // x_init
    if (mxGetField( matlab_model, 0, "x_init" )!=NULL)
        {
        set_x_init = true;
        x_init = mxGetPr( mxGetField( matlab_model, 0, "x_init" ) );
        }




    // ext_qp_res
    if (set_nlp_solver_ext_qp_res)
        {
        ocp_nlp_opts_set(config, opts, "ext_qp_res", &nlp_solver_ext_qp_res);
        }
    // qp_solver_iter_max TODO only for hpipm !!!
    if (set_qp_solver_iter_max)
        {
        ocp_nlp_opts_set(config, opts, "qp_iter_max", &qp_solver_iter_max);
        }
    // qp_solver_cond_N
    if (set_qp_solver_cond_N)
        {
        ocp_nlp_opts_set(config, opts, "qp_cond_N", &qp_solver_cond_N);
        }
    // qp_solver_cond_ric alg
    if (set_qp_solver_cond_N)
        {
        ocp_nlp_opts_set(config, opts, "qp_cond_ric_alg", &qp_solver_cond_ric_alg);
        }
    // qp_solver_ric_alg TODO only for hpipm !!!
    if (set_qp_solver_ric_alg)
        {
        ocp_nlp_opts_set(config, opts, "qp_ric_alg", &qp_solver_ric_alg);
        }
    // qp_solver_warm_start
    if (set_qp_solver_warm_start)
        {
        ocp_nlp_opts_set(config, opts, "qp_warm_start", &qp_solver_warm_start);
        }
    // sim_method_num_stages
    if (set_sim_method_num_stages)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_dynamics_opts_set(config, opts, ii, "num_stages", &sim_method_num_stages);
            }
        }
    // sim_method_num_steps
    if (set_sim_method_num_steps)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_dynamics_opts_set(config, opts, ii, "num_steps", &sim_method_num_steps);
            }
        }




    /* in */

    ocp_nlp_in *in = ocp_nlp_in_create(config, dims);

    // param_scheme_shooting_nodes
    double *param_scheme_shooting_nodes;    bool set_param_scheme_shooting_nodes = false;
    if (mxGetField( matlab_opts, 0, "param_scheme_shooting_nodes" )!=NULL)
	{
        set_param_scheme_shooting_nodes = true;
        param_scheme_shooting_nodes = mxGetPr( mxGetField( matlab_opts, 0, "param_scheme_shooting_nodes" ) );
	}

    // shooting nodes
    // parametrization scheme
    char *param_scheme = mxArrayToString( mxGetField( matlab_opts, 0, "param_scheme" ) );
    if (!strcmp(param_scheme, "multiple_shooting_unif_grid"))
    {
        double Ts = T/N;
        for (ii=0; ii<N; ii++)
       {
            ocp_nlp_in_set(config, dims, in, ii, "Ts", &Ts);
            ocp_nlp_cost_model_set(config, dims, in, ii, "scaling", &Ts);
       }
    }
    else if (!strcmp(param_scheme, "multiple_shooting"))
    {
        if (!set_param_scheme_shooting_nodes)
        {
            mexPrintf("\nerror: ocp_create: param_scheme_shooting_nodes not set for param_scheme multiple_shooting!\n");
            return;
        }
        double scale = T/(param_scheme_shooting_nodes[N]-param_scheme_shooting_nodes[0]);
        for (ii=0; ii<N; ii++)
       {
            double Ts = scale*(param_scheme_shooting_nodes[ii+1]-param_scheme_shooting_nodes[ii]);
            ocp_nlp_in_set(config, dims, in, ii, "Ts", &Ts);
            ocp_nlp_cost_model_set(config, dims, in, ii, "scaling", &Ts);
       }
    }
    else
    {
        mexPrintf("\nerror: ocp_create: param_scheme not supported: %s\n");
        return;
    }

    // cost: ls
    if (!strcmp(cost_type, "linear_ls"))
        {
        // lagrange term
        if (set_Vu)
            {
            for (ii=0; ii<N; ii++)
                {
                ocp_nlp_cost_model_set(config, dims, in, ii, "Vu", Vu);
                }
            }
        if (set_Vx)
            {
            for (ii=0; ii<N; ii++)
                {
                ocp_nlp_cost_model_set(config, dims, in, ii, "Vx", Vx);
                }
            }
        if (set_W)
            {
            for (ii=0; ii<N; ii++)
                {
                ocp_nlp_cost_model_set(config, dims, in, ii, "W", W);
                }
            }
        if (set_yr)
            {
            for (ii=0; ii<N; ii++)
                {
                ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yr);
                }
            }
        }
    if (!strcmp(cost_type_e, "linear_ls"))
        {
        // mayer term
        if (set_Vx_e)
            {
            ocp_nlp_cost_model_set(config, dims, in, N, "Vx", Vx_e);
            }
        if (set_W_e)
            {
            ocp_nlp_cost_model_set(config, dims, in, N, "W", W_e);
            }
        if (set_yr_e)
            {
            ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yr_e);
            }
        }
    // cost: nls
    if (!strcmp(cost_type, "nonlinear_ls"))
        {
        // lagrange term
        if (set_W)
            {
            for (ii=0; ii<N; ii++)
                {
                ocp_nlp_cost_model_set(config, dims, in, ii, "W", W);
                }
            }
        if (set_yr)
            {
            for (ii=0; ii<N; ii++)
                {
                ocp_nlp_cost_model_set(config, dims, in, ii, "y_ref", yr);
                }
            }
        }
    if (!strcmp(cost_type_e, "nonlinear_ls"))
        {
        // mayer term
        if (set_W_e)
            {
            ocp_nlp_cost_model_set(config, dims, in, N, "W", W_e);
            }
        if (set_yr_e)
            {
            ocp_nlp_cost_model_set(config, dims, in, N, "y_ref", yr_e);
            }
        }
    // slacks
    if (set_Z)
        {
        d_ptr = malloc(ns*sizeof(double));
        for (ii=0; ii<ns; ii++)
            {
            d_ptr[ii] = Z[ii+ns*ii];
            }
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_cost_model_set(config, dims, in, ii, "Z", d_ptr);
            }
        free(d_ptr);
        }
    if (set_Z_e)
        {
        d_ptr = malloc(ns*sizeof(double));
        for (ii=0; ii<ns; ii++)
            {
            d_ptr[ii] = Z_e[ii+ns*ii];
            }
        ocp_nlp_cost_model_set(config, dims, in, N, "Z", d_ptr);
        free(d_ptr);
        }
    if (set_Zl)
        {
        d_ptr = malloc(ns*sizeof(double));
        for (ii=0; ii<ns; ii++)
            {
            d_ptr[ii] = Zl[ii+ns*ii];
            }
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_cost_model_set(config, dims, in, ii, "Zl", d_ptr);
            }
        free(d_ptr);
        }
    if (set_Zl_e)
        {
        d_ptr = malloc(ns*sizeof(double));
        for (ii=0; ii<ns; ii++)
            {
            d_ptr[ii] = Zl_e[ii+ns*ii];
            }
        ocp_nlp_cost_model_set(config, dims, in, N, "Zl", d_ptr);
        free(d_ptr);
        }
    if (set_Zu)
        {
        d_ptr = malloc(ns*sizeof(double));
        for (ii=0; ii<ns; ii++)
            {
            d_ptr[ii] = Zu[ii+ns*ii];
            }
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_cost_model_set(config, dims, in, ii, "Zu", d_ptr);
            }
        free(d_ptr);
        }
    if (set_Zu_e)
        {
        d_ptr = malloc(ns*sizeof(double));
        for (ii=0; ii<ns; ii++)
            {
            d_ptr[ii] = Zu_e[ii+ns*ii];
            }
        ocp_nlp_cost_model_set(config, dims, in, N, "Zu", d_ptr);
        free(d_ptr);
        }
    if (set_z)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_cost_model_set(config, dims, in, ii, "z", z);
            }
        }
    if (set_z_e)
        {
        ocp_nlp_cost_model_set(config, dims, in, N, "z", z_e);
        }
    if (set_zl)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_cost_model_set(config, dims, in, ii, "zl", zl);
            }
        }
    if (set_zl_e)
        {
        ocp_nlp_cost_model_set(config, dims, in, N, "zl", zl_e);
        }
    if (set_zu)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_cost_model_set(config, dims, in, ii, "zu", zu);
            }
        }
    if (set_zu_e)
        {
        ocp_nlp_cost_model_set(config, dims, in, N, "zu", zu_e);
        }

    // constraints: bgh

    double acados_inf = 1e8;

    // x0 is always bounded on all components !!!
    i_ptr = malloc(nx*sizeof(int));
    for (ii=0; ii<nx; ii++)
        {
        i_ptr[ii] = ii;
        }
    ocp_nlp_constraints_model_set(config, dims, in, 0, "idxbx", i_ptr);
    free(i_ptr);
    if (set_x0)
        {
        ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", x0);
        ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", x0);
        }
    else
        {
        d_ptr = malloc(nx*sizeof(double));
        for (ii=0; ii<nx; ii++)
            {
            d_ptr[ii] = - acados_inf;
            }
        ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", d_ptr);
        for (ii=0; ii<nx; ii++)
            {
            d_ptr[ii] = acados_inf;
            }
        ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", d_ptr);
        free(d_ptr);
        }
    tmp_idx = malloc(nbx*sizeof(int));
    if (set_Jbx)
        {
        for (ii=0; ii<nbx; ii++)
            {
            idx = -1;
            for (jj=0; jj<nx; jj++)
                {
                if (Jbx[ii+nbx*jj]!=0.0)
                    {
                    tmp_idx[ii] = jj;
                    idx = jj;
                    }
                }
            }
        ii = 1;
        for (; ii<=N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbx", tmp_idx);
            }
        }
    if (set_lbx)
        {
        if (!set_x0)
            {
            d_ptr = malloc(nx*sizeof(double));
            for (ii=0; ii<nx; ii++)
                {
                d_ptr[ii] = - acados_inf;
                }
            for (ii=0; ii<nbx; ii++)
                {
                d_ptr[tmp_idx[ii]] = lbx[ii];
                }
            ocp_nlp_constraints_model_set(config, dims, in, 0, "lbx", d_ptr);
            free(d_ptr);
            }
        ii = 1;
        for (; ii<=N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lbx", lbx);
            }
        }
    if (set_ubx)
        {
        if (!set_x0)
            {
            d_ptr = malloc(nx*sizeof(double));
            for (ii=0; ii<nx; ii++)
                {
                d_ptr[ii] = acados_inf;
                }
            for (ii=0; ii<nbx; ii++)
                {
                d_ptr[tmp_idx[ii]] = ubx[ii];
                }
            ocp_nlp_constraints_model_set(config, dims, in, 0, "ubx", d_ptr);
            free(d_ptr);
            }
        ii = 1;
        for (; ii<=N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "ubx", ubx);
            }
        }
    free(tmp_idx);
    if (set_Jbu)
        {
        i_ptr = malloc(nbu*sizeof(int));
        for (ii=0; ii<nbu; ii++)
            {
            idx = -1;
            for (jj=0; jj<nu; jj++)
                {
                if (Jbu[ii+nbu*jj]!=0.0)
                    {
                    i_ptr[ii] = jj;
                    idx = jj;
                    }
                }
            }
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxbu", i_ptr);
            }
        free(i_ptr);
        }
    if (set_lbu)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lbu", lbu);
            }
        }
    if (set_ubu)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "ubu", ubu);
            }
        }
    if (set_C)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "C", C);
            }
        }
    if (set_D)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "D", D);
            }
        }
    if (set_lg)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lg", lg);
            }
        }
    if (set_ug)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "ug", ug);
            }
        }
    if (set_C_e)
        {
        ocp_nlp_constraints_model_set(config, dims, in, N, "C", C_e);
        }
    if (set_lg_e)
        {
        ocp_nlp_constraints_model_set(config, dims, in, N, "lg", lg_e);
        }
    if (set_ug_e)
        {
        ocp_nlp_constraints_model_set(config, dims, in, N, "ug", ug_e);
        }
    if (set_lh)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "lh", lh);
            }
        }
    if (set_lh_e)
        {
        ocp_nlp_constraints_model_set(config, dims, in, N, "lh", lh_e);
        }
    if (set_uh)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "uh", uh);
            }
        }
    if (set_uh_e)
        {
        ocp_nlp_constraints_model_set(config, dims, in, N, "uh", uh_e);
        }
    if (set_Jsbu)
        {
        i_ptr = malloc(nsbu*sizeof(int));
        for (ii=0; ii<nsbu; ii++)
            {
            idx = -1;
            for (jj=0; jj<nbu; jj++)
                {
                if (Jsbu[jj+nbu*ii]!=0.0)
                    {
                    i_ptr[ii] = jj;
                    idx = jj;
                    }
                }
            }
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsbu", i_ptr);
            }
        free(i_ptr);
        }
//    if (set_lsbu)
//        {
//        for (ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsbu", lsbu);
//            }
//        }
//    if (set_usbu)
//        {
//        for (ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "usbu", usbu);
//            }
//        }
    if (set_Jsbx)
        {
        i_ptr = malloc(nsbx*sizeof(int));
        for (ii=0; ii<nsbx; ii++)
            {
            idx = -1;
            for (jj=0; jj<nbx; jj++)
                {
                if (Jsbx[jj+nbx*ii]!=0.0)
                    {
                    i_ptr[ii] = jj;
                    idx = jj;
                    }
                }
            }
//        for (ii=0; ii<=N; ii++)
        for (ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsbx", i_ptr);
            }
        free(i_ptr);
        }
//    if (set_lsbx)
//        {
//        for (ii=0; ii<=N; ii++)
//        for (ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsbx", lsbx);
//            }
//        }
//    if (set_usbx)
//        {
//        for (ii=0; ii<=N; ii++)
//        for (ii=1; ii<=N; ii++) // TODO stage 0 !!!!!!!!!!
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "usbx", usbx);
//            }
//        }
    if (set_Jsg)
        {
        i_ptr = malloc(nsg*sizeof(int));
        for (ii=0; ii<nsg; ii++)
            {
            idx = -1;
            for (jj=0; jj<ng; jj++)
                {
                if (Jsg[jj+ng*ii]!=0.0)
                    {
                    i_ptr[ii] = jj;
                    idx = jj;
                    }
                }
            }
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsg", i_ptr);
            }
        free(i_ptr);
        }
//    if (set_lsg)
//        {
//        for (ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsg", lsg);
//            }
//        }
//    if (set_usg)
//        {
//        for (ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "usg", usg);
//            }
//        }
    if (set_Jsg_e)
        {
        i_ptr = malloc(nsg_e*sizeof(int));
        for (ii=0; ii<nsg_e; ii++)
            {
            idx = -1;
            for (jj=0; jj<ng_e; jj++)
                {
                if (Jsg_e[jj+ng_e*ii]!=0.0)
                    {
                    i_ptr[ii] = jj;
                    idx = jj;
                    }
                }
            }
        ocp_nlp_constraints_model_set(config, dims, in, N, "idxsg", i_ptr);
        free(i_ptr);
        }
//    if (set_lsg_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "lsg", lsg_e);
//        }
//    if (set_usg_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "usg", usg_e);
//        }
    if (set_Jsh)
        {
        i_ptr = malloc(nsh*sizeof(int));
        for (ii=0; ii<nsh; ii++)
            {
            idx = -1;
            for (jj=0; jj<nh; jj++)
                {
                if (Jsh[jj+nh*ii]!=0.0)
                    {
                    i_ptr[ii] = jj;
                    idx = jj;
                    }
                }
            }
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_constraints_model_set(config, dims, in, ii, "idxsh", i_ptr);
            }
        free(i_ptr);
        }
//    if (set_lsh)
//        {
//        for (ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "lsh", lsh);
//            }
//        }
//    if (set_ush)
//        {
//        for (ii=0; ii<N; ii++)
//            {
//            ocp_nlp_constraints_model_set(config, dims, in, ii, "ush", ush);
//            }
//        }
    if (set_Jsh_e)
        {
        i_ptr = malloc(nsh_e*sizeof(int));
        for (ii=0; ii<nsh_e; ii++)
            {
            idx = -1;
            for (jj=0; jj<nh_e; jj++)
                {
                if (Jsh_e[jj+nh_e*ii]!=0.0)
                    {
                    i_ptr[ii] = jj;
                    idx = jj;
                    }
                }
            }
        ocp_nlp_constraints_model_set(config, dims, in, N, "idxsh", i_ptr);
        free(i_ptr);
        }
//    if (set_lsh_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "lsh", lsh_e);
//        }
//    if (set_ush_e)
//        {
//        ocp_nlp_constraints_model_set(config, dims, in, N, "ush", ush_e);
//        }



    /* out */

    ocp_nlp_out *out = ocp_nlp_out_create(config, dims);

    if (set_x_init)
        {
        for (ii=0; ii<=N; ii++)
            {
            ocp_nlp_out_set(config, dims, out, ii, "x", x_init+ii*nx);
            }
        }
    else // initialize to zero
        {
        double *x_init = calloc(nx, sizeof(double));
        for (ii=0; ii<=N; ii++)
            {
            ocp_nlp_out_set(config, dims, out, ii, "x", x_init);
            }
        free(x_init);
        }
    if (set_u_init)
        {
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_out_set(config, dims, out, ii, "u", u_init+ii*nu);
            }
        }
    else // initialize to zero
        {
        double *u_init = calloc(nu, sizeof(double));
        for (ii=0; ii<N; ii++)
            {
            ocp_nlp_out_set(config, dims, out, ii, "u", u_init);
            }
        free(u_init);
        }


    /* solver */
    ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, opts);


    /* sens_out */
    ocp_nlp_out *sens_out = ocp_nlp_out_create(config, dims);


    /* populate output struct */
    // plan
    mxArray *plan_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(plan_mat);
    l_ptr[0] = (long long) plan;
    mxSetField(plhs[0], 0, "plan", plan_mat);

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

    // sens_out
    mxArray *sens_out_mat  = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    l_ptr = mxGetData(sens_out_mat);
    l_ptr[0] = (long long) sens_out;
    mxSetField(plhs[0], 0, "sens_out", sens_out_mat);


    return;

}



#define S_FUNCTION_NAME engine_nmpc
#define S_FUNCTION_LEVEL 2

#ifndef MATLAB_MEX_FILE
#include <brtenv.h>
#define printf(...) msg_info_printf(MSG_SM_USER, 0, __VA_ARGS__);
#endif

#define NUM_STAGES 20
#define NUM_STATES 4
#define NUM_CONTROLS 2
#define HORIZON_LENGTH 1.0

#define ACADOS_WITH_QPOASES
#define ACADOS_WITH_HPMPC

#include "simstruc.h"

#include <stdlib.h>

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/math.h"
#include "acados/utils/print.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_nlp_interface.h"

#include "examples/c/engine_model/engine_impl_dae_fun.h"
#include "examples/c/engine_model/engine_impl_dae_fun_jac_x_xdot_z.h"
#include "examples/c/engine_model/engine_impl_dae_jac_x_xdot_u_z.h"
#include "examples/c/engine_model/engine_ls_cost.h"
#include "examples/c/engine_model/engine_ls_cost_N.h"
#include "examples/c/engine_model/reference.c"

#define CASADI_WORK_FUNCTION_CAT(a) a##_work
#define CASADI_SPARSITY_IN_FUNCTION_CAT(a) a##_sparsity_in
#define CASADI_SPARSITY_OUT_FUNCTION_CAT(a) a##_sparsity_out
#define CASADI_N_IN_FUNCTION_CAT(a) a##_n_in
#define CASADI_N_OUT_FUNCTION_CAT(a) a##_n_out

#define CASADI_WORK_FUNCTION(a) CASADI_WORK_FUNCTION_CAT(a)
#define CASADI_SPARSITY_IN_FUNCTION(a) CASADI_SPARSITY_IN_FUNCTION_CAT(a)
#define CASADI_SPARSITY_OUT_FUNCTION(a) CASADI_SPARSITY_OUT_FUNCTION_CAT(a)
#define CASADI_N_IN_FUNCTION(a) CASADI_N_IN_FUNCTION_CAT(a)
#define CASADI_N_OUT_FUNCTION(a) CASADI_N_OUT_FUNCTION_CAT(a)

external_function_casadi impl_dae_fun, impl_dae_fun_jac_x_xdot_z, impl_dae_jac_x_xdot_u_z, nls_cost_residual, nls_cost_N_residual;

static void mdlInitializeSizes(SimStruct *S)
{

    if (!ssSetNumInputPorts(S, 2)) return;
    ssSetInputPortWidth(S, 0, NUM_STATES);  // x0
    ssSetInputPortWidth(S, 1, 1);           // reference

    ssSetInputPortDirectFeedThrough(S, 0, true);
    ssSetInputPortRequiredContiguous(S, 0, true);
    ssSetInputPortDirectFeedThrough(S, 1, true);
    ssSetInputPortRequiredContiguous(S, 1, true);

    ssSetNumPWork(S, 6);

    if (!ssSetNumOutputPorts(S, 4)) return;
    ssSetOutputPortWidth(S, 0, NUM_CONTROLS);
    ssSetOutputPortWidth(S, 1, NUM_STATES);
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortWidth(S, 3, 1);

    ssSetNumSampleTimes(S, 1);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    int i, j;

    int nx[NUM_STAGES+1], nu[NUM_STAGES+1], nz[NUM_STAGES+1], ny[NUM_STAGES+1], nb[NUM_STAGES+1], nbx[NUM_STAGES+1], nbu[NUM_STAGES+1], ng[NUM_STAGES+1], nh[NUM_STAGES+1], nq[NUM_STAGES+1], ns[NUM_STAGES+1];
    for (i = 0; i <= NUM_STAGES; ++i)
    {
        nx[i] = NUM_STATES;
        nu[i] = NUM_CONTROLS;
        nz[i] = 2;
        ny[i] = 1 + nx[i] + nu[i];
        nbx[i] = nx[i];
        nbu[i] = nu[i];
        nb[i] = nbx[i] + nbu[i];
        ng[i] = 0;
        nh[i] = 0;
        nq[i] = 0;
        ns[i] = 0;
    }

    nu[NUM_STAGES] = 0;
    nz[NUM_STAGES] = 0;
    ny[NUM_STAGES] = 1 + nx[NUM_STAGES];
    nbu[NUM_STAGES] = 0;
    nb[NUM_STAGES] = nbx[NUM_STAGES] + nbu[NUM_STAGES];

    int idxb[nbu[0]+nbx[0]];
    for (i = 0; i < nbu[0]+nbx[0]; ++i)
        idxb[i] = i;

    // sampling time (s)
    double T = 0.05;

    // x: u1, u2, xD1, xD2
    double x0[] = {50, 50, 1.14275, 1.53787};
    // z: xA1, xA2
    static double z0[] = {1.28976, 1.78264};
    static struct blasfeo_dvec z0_dvec;
    blasfeo_allocate_dvec(2, &z0_dvec);
    blasfeo_pack_dvec(2, z0, &z0_dvec, 0);
    // u: u1_r, u2_r
    double u[] = {0, 0};

    // reference
    double y_ref[] = {1.500, 50, 50, 1.14275, 1.53787, 0, 0};

    // weighting matrices
    double W[ny[0]*ny[0]];
    for (i = 0; i < ny[0]*ny[0]; ++i)
        W[i] = 0;
    W[0*(ny[0]+1)] = 1000;
    W[1*(ny[0]+1)] = 1e-3;
    W[2*(ny[0]+1)] = 1e-3;
    W[3*(ny[0]+1)] = 1e-3;
    W[4*(ny[0]+1)] = 1e-3;
    W[5*(ny[0]+1)] = 0.1e-3;
    W[6*(ny[0]+1)] = 0.1e-3;

    double W_N[ny[NUM_STAGES]*ny[NUM_STAGES]];
    for (i = 0; i < ny[NUM_STAGES]*ny[NUM_STAGES]; ++i)
        W_N[i] = 0;
    W_N[0*(ny[NUM_STAGES]+1)] = 1000;
    W_N[1*(ny[NUM_STAGES]+1)] = 1e-3;
    W_N[2*(ny[NUM_STAGES]+1)] = 1e-3;
    W_N[3*(ny[NUM_STAGES]+1)] = 1e-3;
    W_N[4*(ny[NUM_STAGES]+1)] = 1e-3;

    double lb_0[] = {-10000, -10000, 50, 50, 1.14275, 1.53787};
    double ub_0[] = {+10000, +10000, 50, 50, 1.14275, 1.53787};

    double lb[] = {-10000, -10000, 0, 0, 0.5, 0.5};
    double ub[] = {+10000, +10000, 100, 100, 1.757, 2.125};

    double lb_N[] = {0, 0, 0.5, 0.5};
    double ub_N[] = {100, 100, 1.757, 2.125};

    ocp_nlp_plan *plan = ocp_nlp_plan_create(NUM_STAGES);

	plan->nlp_solver = SQP_GN;

	for (i = 0; i <= NUM_STAGES; i++)
		plan->nlp_cost[i] = NONLINEAR_LS;

	plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;

	for (i = 0; i < NUM_STAGES; i++)
    {
		plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
		plan->sim_solver_plan[i].sim_solver = IRK;
	}

	for (i = 0; i <= NUM_STAGES; i++)
		plan->nlp_constraints[i] = BGH;

	ocp_nlp_config *config = ocp_nlp_config_create(*plan);

    // implicit dae
    impl_dae_fun.casadi_fun = &engine_impl_dae_fun;
    impl_dae_fun.casadi_work = &engine_impl_dae_fun_work;
    impl_dae_fun.casadi_sparsity_in = &engine_impl_dae_fun_sparsity_in;
    impl_dae_fun.casadi_sparsity_out = &engine_impl_dae_fun_sparsity_out;
    impl_dae_fun.casadi_n_in = &engine_impl_dae_fun_n_in;
    impl_dae_fun.casadi_n_out = &engine_impl_dae_fun_n_out;
    external_function_casadi_create(&impl_dae_fun);

    impl_dae_fun_jac_x_xdot_z.casadi_fun = &engine_impl_dae_fun_jac_x_xdot_z;
    impl_dae_fun_jac_x_xdot_z.casadi_work = &engine_impl_dae_fun_jac_x_xdot_z_work;
    impl_dae_fun_jac_x_xdot_z.casadi_sparsity_in = &engine_impl_dae_fun_jac_x_xdot_z_sparsity_in;
    impl_dae_fun_jac_x_xdot_z.casadi_sparsity_out = &engine_impl_dae_fun_jac_x_xdot_z_sparsity_out;
    impl_dae_fun_jac_x_xdot_z.casadi_n_in = &engine_impl_dae_fun_jac_x_xdot_z_n_in;
    impl_dae_fun_jac_x_xdot_z.casadi_n_out = &engine_impl_dae_fun_jac_x_xdot_z_n_out;
    external_function_casadi_create(&impl_dae_fun_jac_x_xdot_z);

    impl_dae_jac_x_xdot_u_z.casadi_fun = &engine_impl_dae_jac_x_xdot_u_z;
    impl_dae_jac_x_xdot_u_z.casadi_work = &engine_impl_dae_jac_x_xdot_u_z_work;
    impl_dae_jac_x_xdot_u_z.casadi_sparsity_in = &engine_impl_dae_jac_x_xdot_u_z_sparsity_in;
    impl_dae_jac_x_xdot_u_z.casadi_sparsity_out = &engine_impl_dae_jac_x_xdot_u_z_sparsity_out;
    impl_dae_jac_x_xdot_u_z.casadi_n_in = &engine_impl_dae_jac_x_xdot_u_z_n_in;
    impl_dae_jac_x_xdot_u_z.casadi_n_out = &engine_impl_dae_jac_x_xdot_u_z_n_out;
    external_function_casadi_create(&impl_dae_jac_x_xdot_u_z);

    nls_cost_residual.casadi_fun = &engine_ls_cost;
    nls_cost_residual.casadi_work = &engine_ls_cost_work;
    nls_cost_residual.casadi_sparsity_in = &engine_ls_cost_sparsity_in;
    nls_cost_residual.casadi_sparsity_out = &engine_ls_cost_sparsity_out;
    nls_cost_residual.casadi_n_in = &engine_ls_cost_n_in;
    nls_cost_residual.casadi_n_out = &engine_ls_cost_n_out;
    external_function_casadi_create(&nls_cost_residual);

    nls_cost_N_residual.casadi_fun = &engine_ls_cost_N;
    nls_cost_N_residual.casadi_work = &engine_ls_cost_N_work;
    nls_cost_N_residual.casadi_sparsity_in = &engine_ls_cost_N_sparsity_in;
    nls_cost_N_residual.casadi_sparsity_out = &engine_ls_cost_N_sparsity_out;
    nls_cost_N_residual.casadi_n_in = &engine_ls_cost_N_n_in;
    nls_cost_N_residual.casadi_n_out = &engine_ls_cost_N_n_out;
    external_function_casadi_create(&nls_cost_N_residual);

    // dimensions
    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);

    ocp_nlp_dims_set_opt_vars(config, dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(config, dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", ns);

	for (int i = 0; i <= N; i++)
    {
        ocp_nlp_dims_set_cost(config, dims, i, "ny", &ny[i]);

        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", &nh[i]);
    }

    // in
	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);
	for (i = 0; i < NUM_STAGES; ++i)
    	nlp_in->Ts[i] = T;

    // cost
    ocp_nlp_cost_nls_model **cost = (ocp_nlp_cost_nls_model **) nlp_in->cost;

	for (i = 0; i < NUM_STAGES; ++i) {
        cost[i]->nls_res_jac = (external_function_generic *) &nls_cost_residual;
        cost[i]->nls_hess = NULL;
        blasfeo_pack_dvec(ny[i], y_ref, &cost[i]->y_ref, 0);
        blasfeo_pack_dmat(ny[i], ny[i], W, ny[i], &cost[i]->W, 0, 0);
    }

    cost[NUM_STAGES]->nls_res_jac = (external_function_generic *) &nls_cost_N_residual;
    cost[NUM_STAGES]->nls_hess = NULL;
    blasfeo_pack_dvec(ny[NUM_STAGES], y_ref, &cost[NUM_STAGES]->y_ref, 0);
    blasfeo_pack_dmat(ny[NUM_STAGES], ny[NUM_STAGES], W_N, ny[NUM_STAGES], &cost[NUM_STAGES]->W, 0, 0);

    // dynamics
    for (i = 0; i < NUM_STAGES; ++i)
    {
        if(ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun", &impl_dae_fun)) exit(1);
        if(ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_dae_fun_jac_x_xdot_z)) exit(1);
        if(ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_dae_jac_x_xdot_u_z)) exit(1);
    }

    // bounds
	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;
	ocp_nlp_constraints_bgh_dims **constraints_dims = (ocp_nlp_constraints_bgh_dims **) dims->constraints;

    ocp_nlp_constraints_bounds_set(config, dims, nlp_in, 0, "lb", lb_0);
    ocp_nlp_constraints_bounds_set(config, dims, nlp_in, 0, "ub", ub_0);
    for (i = 0; i < nb[0]; ++i)
        constraints[0]->idxb[i] = idxb[i];

    for (i = 1; i < NUM_STAGES; ++i)
    {
        ocp_nlp_constraints_bounds_set(config, dims, nlp_in, i, "lb", lb);
        ocp_nlp_constraints_bounds_set(config, dims, nlp_in, i, "ub", ub);
        for (j = 0; j < nb[i]; ++j)
            constraints[i]->idxb[j] = idxb[j];
    }

    ocp_nlp_constraints_bounds_set(config, dims, nlp_in, NUM_STAGES, "lb", lb_N);
    ocp_nlp_constraints_bounds_set(config, dims, nlp_in, NUM_STAGES, "ub", ub_N);

    for (i = 0; i < nb[NUM_STAGES]; ++i)
        constraints[NUM_STAGES]->idxb[i] = idxb[i];

    // options
    ocp_nlp_sqp_opts *nlp_opts = ocp_nlp_opts_create(config, dims);
    nlp_opts->maxIter = 1;

    // out
    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

    // solver
	ocp_nlp_solver *solver = ocp_nlp_create(config, dims, nlp_opts);

    // initialize
    for (i = 0; i < NUM_STAGES; ++i)
    {
        blasfeo_pack_dvec(nu[i], u, nlp_out->ux+i, 0);
        blasfeo_pack_dvec(nx[i], x0, nlp_out->ux+i, nu[i]);
        nlp_out->z = &z0_dvec;
    }
    blasfeo_pack_dvec(nx[NUM_STAGES], x0, nlp_out->ux+NUM_STAGES, nu[NUM_STAGES]);

    ssGetPWork(S)[0] = (void *) dims;
    ssGetPWork(S)[1] = (void *) nlp_in;
    ssGetPWork(S)[2] = (void *) nlp_out;
    ssGetPWork(S)[3] = (void *) nlp_opts;
    ssGetPWork(S)[4] = (void *) solver;
    ssGetPWork(S)[5] = (void *) config;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    int j;

    ocp_nlp_dims *nlp_dims = (ocp_nlp_dims *) ssGetPWork(S)[0];
    ocp_nlp_in *nlp_in = (ocp_nlp_in *) ssGetPWork(S)[1];
    ocp_nlp_out *nlp_out = (ocp_nlp_out *) ssGetPWork(S)[2];
    ocp_nlp_solver *nlp_solver = (ocp_nlp_solver *) ssGetPWork(S)[4];

    ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;
	ocp_nlp_constraints_bgh_dims **constraints_dims = (ocp_nlp_constraints_bgh_dims **) nlp_dims->constraints;
    ocp_nlp_cost_nls_model **cost = (ocp_nlp_cost_nls_model **) nlp_in->cost;

    const double *x0 = ssGetInputPortRealSignal(S, 0);
    const double *reference = ssGetInputPortRealSignal(S, 1);

    // bounds
    double lb_0[] = {-10000, -10000, 50, 50, 1.14275, 1.53787};
    double ub_0[] = {+10000, +10000, 50, 50, 1.14275, 1.53787};
    for (j = 0; j < NUM_STATES; ++j) {
        lb_0[NUM_CONTROLS+j] = x0[j];
        ub_0[NUM_CONTROLS+j] = x0[j];
    } 

    ocp_nlp_constraints_bounds_set(config, dims, nlp_in, 0, "lb", lb_0);
    ocp_nlp_constraints_bounds_set(config, dims, nlp_in, 0, "ub", ub_0);
    
    for (j = 0; j <= NUM_STAGES; ++j)
        BLASFEO_DVECEL(&cost[j]->y_ref, 0) = *reference;

    int status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);
    
    double *u0_opt = ssGetOutputPortRealSignal(S, 0);
    double *x1 = ssGetOutputPortRealSignal(S, 1);
    double *status_out = ssGetOutputPortRealSignal(S, 2);
    double *comp_time = ssGetOutputPortRealSignal(S, 3);

    blasfeo_unpack_dvec(NUM_CONTROLS, &nlp_out->ux[0], 0, u0_opt);
    blasfeo_unpack_dvec(NUM_STATES, &nlp_out->ux[1], NUM_CONTROLS, x1);
    *status_out = (double) status;
    *comp_time = nlp_out->total_time;
}

static void mdlTerminate(SimStruct *S)
{
    free(ssGetPWork(S)[0]);
    free(ssGetPWork(S)[1]);
    free(ssGetPWork(S)[2]);
    free(ssGetPWork(S)[3]);
    free(ssGetPWork(S)[4]);
    free(ssGetPWork(S)[5]);
}

#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c"  /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif

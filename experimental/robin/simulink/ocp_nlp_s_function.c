
#define S_FUNCTION_NAME ocp_nlp_s_function
#define S_FUNCTION_LEVEL 2

#ifndef MATLAB_MEX_FILE
#include <brtenv.h>
#define printf(...) msg_info_printf(MSG_SM_USER, 0, __VA_ARGS__);
#endif

// Fill in problem specific properties >>>>>
#define NUM_STAGES 5
#define HORIZON_LENGTH 1.0
#define MODEL_NAME vdeFun
#include "../examples/c/crane_model/crane_model.h"
// <<<<< until here

#include "simstruc.h"

#include "blasfeo/include/blasfeo_d_aux.h"

#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/utils/print.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/options_interface.h"

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

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 2);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
    {
        return; /* Parameter mismatch will be reported by Simulink */
    }
    const mxArray *Q = ssGetSFcnParam(S, 0);
    const mxArray *R = ssGetSFcnParam(S, 1);

    int nx = mxGetM(Q);
    int nu = mxGetM(R);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, nx);  // x0

    ssSetInputPortDirectFeedThrough(S, 0, true);
    ssSetInputPortRequiredContiguous(S, 0, true);

    ssSetNumPWork(S, 5);

    if (!ssSetNumOutputPorts(S, 3)) return;
    ssSetOutputPortWidth(S, 0, nu);
    ssSetOutputPortWidth(S, 1, nx);
    ssSetOutputPortWidth(S, 2, 1);

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
    const mxArray *Q = ssGetSFcnParam(S, 0);
    const mxArray *R = ssGetSFcnParam(S, 1);

    int num_states = mxGetM(Q);
    ssSetInputPortWidth(S, 0, num_states);
    int num_controls = mxGetM(R);
    ssSetOutputPortWidth(S, 0, num_controls);

    int nx[NUM_STAGES + 1], nu[NUM_STAGES + 1], ny[NUM_STAGES + 1], nb[NUM_STAGES + 1],
        nbx[NUM_STAGES + 1], nbu[NUM_STAGES + 1], ng[NUM_STAGES + 1], nh[NUM_STAGES + 1],
        ns[NUM_STAGES + 1], nq[NUM_STAGES + 1];

    for (int i = 0; i < NUM_STAGES; ++i)
    {
        nx[i] = num_states;
        nu[i] = num_controls;
        ny[i] = num_states + num_controls;
        nbx[i] = 0;
        nbu[i] = 0;
        nb[i] = nbx[i] + nbu[i];
        ng[i] = 0;
        nh[i] = 0;
        ns[i] = 0;
        nq[i] = 0;
    }

    nbx[0] = num_states;
    nb[0] = nbx[0] + nbu[0];

    nx[NUM_STAGES] = num_states;
    nu[NUM_STAGES] = 0;
    ny[NUM_STAGES] = num_states;
    nbx[NUM_STAGES] = 0;
    nbu[NUM_STAGES] = 0;
    ng[NUM_STAGES] = 0;
    nh[NUM_STAGES] = 0;
    ns[NUM_STAGES] = 0;
    nq[NUM_STAGES] = 0;

    ocp_nlp_solver_plan *plan = ocp_nlp_plan_create(NUM_STAGES);
    plan->nlp_solver = SQP;

    for (int i = 0; i < NUM_STAGES; i++)
    {
        plan->nlp_cost[i] = LINEAR_LS;
        plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        plan->sim_solver_plan[i].sim_solver = ERK;
    }
    plan->nlp_cost[NUM_STAGES] = LINEAR_LS;

    plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;

    ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan);

    ocp_nlp_dims *nlp_dims = ocp_nlp_dims_create(config);
    ocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, ns, nq, nlp_dims);

    external_function_casadi *expl_vde_for = malloc(NUM_STAGES * sizeof(external_function_casadi));
    for (int i = 0; i < NUM_STAGES; ++i)
    {
        expl_vde_for[i].casadi_fun = &MODEL_NAME;
        expl_vde_for[i].casadi_work = &CASADI_WORK_FUNCTION(MODEL_NAME);
        expl_vde_for[i].casadi_sparsity_in = &CASADI_SPARSITY_IN_FUNCTION(MODEL_NAME);
        expl_vde_for[i].casadi_sparsity_out = &CASADI_SPARSITY_OUT_FUNCTION(MODEL_NAME);
        expl_vde_for[i].casadi_n_in = &CASADI_N_IN_FUNCTION(MODEL_NAME);
        expl_vde_for[i].casadi_n_out = &CASADI_N_OUT_FUNCTION(MODEL_NAME);
    }

    external_function_casadi_create_array(NUM_STAGES, expl_vde_for);

    ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, nlp_dims);

    for (int i = 0; i < NUM_STAGES; ++i) nlp_in->Ts[i] = HORIZON_LENGTH / NUM_STAGES;

    ocp_nlp_cost_ls_model *stage_cost_ls;
    for (int i = 0; i <= NUM_STAGES; ++i)
    {
        stage_cost_ls = (ocp_nlp_cost_ls_model *) nlp_in->cost[i];
        // Cyt
        blasfeo_dgese(nu[i] + nx[i], ny[i], 0.0, &stage_cost_ls->Cyt, 0, 0);
        for (int j = 0; j < nu[i]; j++) BLASFEO_DMATEL(&stage_cost_ls->Cyt, j, nx[i] + j) = 1.0;
        for (int j = 0; j < nx[i]; j++) BLASFEO_DMATEL(&stage_cost_ls->Cyt, nu[i] + j, j) = 1.0;

        // W
        blasfeo_dgese(ny[i], ny[i], 0.0, &stage_cost_ls->W, 0, 0);
        blasfeo_pack_dmat(nx[i], nx[i], (double *) mxGetData(Q), nx[i], &stage_cost_ls->W, 0, 0);
        blasfeo_pack_dmat(nu[i], nu[i], (double *) mxGetData(R), nu[i], &stage_cost_ls->W, nx[i],
                          nx[i]);

        // y_ref
        blasfeo_dvecse(nx[i], 0.0, &stage_cost_ls->y_ref, 0);
        blasfeo_dvecse(nu[i], 0.0, &stage_cost_ls->y_ref, nx[i]);
    }

    for (int i = 0; i < NUM_STAGES; ++i)
    {
        ocp_nlp_dynamics_model_set(config, nlp_in, i, "expl_vde_for", &expl_vde_for[i]);
    }

    ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_in->constraints;

    for (int i = 0; i < num_states; ++i) constraints[0]->idxb[i] = num_controls + i;

    void *nlp_opts = ocp_nlp_opts_create(config, nlp_dims);

    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, nlp_dims);

    ocp_nlp_solver *nlp_solver = ocp_nlp_create(config, nlp_dims, nlp_opts);

    ssGetPWork(S)[0] = (void *) nlp_dims;
    ssGetPWork(S)[1] = (void *) nlp_in;
    ssGetPWork(S)[2] = (void *) nlp_out;
    ssGetPWork(S)[3] = (void *) nlp_opts;
    ssGetPWork(S)[4] = (void *) nlp_solver;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    ocp_nlp_in *nlp_in = (ocp_nlp_in *) ssGetPWork(S)[1];
    ocp_nlp_out *nlp_out = (ocp_nlp_out *) ssGetPWork(S)[2];
    ocp_nlp_solver *nlp_solver = (ocp_nlp_solver *) ssGetPWork(S)[4];

    const mxArray *Q = ssGetSFcnParam(S, 0);
    const mxArray *R = ssGetSFcnParam(S, 1);

    int nx = mxGetM(Q);
    int nu = mxGetM(R);

    ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_in->constraints;

    const double *x0 = ssGetInputPortRealSignal(S, 0);
    blasfeo_pack_dvec(nx, (double *) x0, &constraints[0]->d, 0);
    blasfeo_pack_dvec(nx, (double *) x0, &constraints[0]->d, nx);

    int status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);

    double *u0_opt = ssGetOutputPortRealSignal(S, 0);
    double *x1 = ssGetOutputPortRealSignal(S, 1);
    double *status_out = ssGetOutputPortRealSignal(S, 2);

    blasfeo_unpack_dvec(nu, &nlp_out->ux[0], 0, u0_opt);
    blasfeo_unpack_dvec(nx, &nlp_out->ux[1], nu, x1);
    *status_out = (double) status;
}

static void mdlTerminate(SimStruct *S)
{
    free(ssGetPWork(S)[0]);
    free(ssGetPWork(S)[1]);
    free(ssGetPWork(S)[2]);
    free(ssGetPWork(S)[3]);
    free(ssGetPWork(S)[4]);
}

#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c"  /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif

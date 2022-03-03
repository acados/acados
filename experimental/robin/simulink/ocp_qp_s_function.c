
#define S_FUNCTION_NAME ocp_qp_s_function
#define S_FUNCTION_LEVEL 2

#ifndef MATLAB_MEX_FILE
#include <brtenv.h>
#define printf(...) msg_info_printf(MSG_SM_USER, 0, __VA_ARGS__);
#endif

#define AC_HORIZON_LENGTH 5

#include "simstruc.h"

#include "acados/utils/print.h"
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/options_interface.h"

static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 4);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
    {
        return; /* Parameter mismatch will be reported by Simulink */
    }
    const mxArray *B = ssGetSFcnParam(S, 3);

    int nx = mxGetM(B);
    int nu = mxGetN(B);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, nx);
    ssSetInputPortDirectFeedThrough(S, 0, true);
    ssSetInputPortRequiredContiguous(S, 0, true);

    ssSetNumPWork(S, 5);

    if (!ssSetNumOutputPorts(S, 3)) return;
    ssSetOutputPortWidth(S, 0, nu);
    ssSetOutputPortWidth(S, 1, 1);
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
    printf("mdlStart: Read parameters");
    const mxArray *Q = ssGetSFcnParam(S, 0);
    const mxArray *R = ssGetSFcnParam(S, 1);
    const mxArray *A = ssGetSFcnParam(S, 2);
    const mxArray *B = ssGetSFcnParam(S, 3);

    int nx = mxGetM(B);
    ssSetInputPortWidth(S, 0, nx);
    int nu = mxGetN(B);
    ssSetOutputPortWidth(S, 0, nu);

    ocp_qp_dims *qp_dims = ocp_qp_dims_create(AC_HORIZON_LENGTH);
    qp_dims->nbx[0] = nx;
    qp_dims->nb[0] = nx;
    for (int i = 0; i < AC_HORIZON_LENGTH; ++i)
    {
        qp_dims->nx[i] = nx;
        qp_dims->nu[i] = nu;
    }
    qp_dims->nx[AC_HORIZON_LENGTH] = nx;

    ocp_qp_solver_plan_t plan = {PARTIAL_CONDENSING_HPIPM};
    ocp_qp_xcond_solver_config *config = ocp_qp_config_create(plan);

    ocp_qp_in *qp_in = ocp_qp_in_create(config, qp_dims);
    for (int i = 0; i < nx; ++i) qp_in->idxb[0][i] = nu + i;

    for (int i = 0; i < AC_HORIZON_LENGTH; ++i)
    {
        d_cvt_colmaj_to_ocp_qp_Q(i, (double *) mxGetData(Q), qp_in);
        d_cvt_colmaj_to_ocp_qp_R(i, (double *) mxGetData(R), qp_in);
        d_cvt_colmaj_to_ocp_qp_A(i, (double *) mxGetData(A), qp_in);
        d_cvt_colmaj_to_ocp_qp_B(i, (double *) mxGetData(B), qp_in);
    }
    d_cvt_colmaj_to_ocp_qp_Q(AC_HORIZON_LENGTH, (double *) mxGetData(Q), qp_in);
    d_cvt_colmaj_to_ocp_qp_R(AC_HORIZON_LENGTH, (double *) mxGetData(R), qp_in);

    ocp_qp_out *qp_out = ocp_qp_out_create(config, qp_dims);

    void *qp_opts = ocp_qp_opts_create(config, qp_dims);

    printf("mdlStart: Create acados solver");
    ocp_qp_solver *qp_solver = ocp_qp_create(config, qp_dims, qp_opts);

    ssGetPWork(S)[0] = (void *) qp_dims;
    ssGetPWork(S)[1] = (void *) qp_in;
    ssGetPWork(S)[2] = (void *) qp_out;
    ssGetPWork(S)[3] = (void *) qp_opts;
    ssGetPWork(S)[4] = (void *) qp_solver;
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    ocp_qp_in *qp_in = (ocp_qp_in *) ssGetPWork(S)[1];
    ocp_qp_out *qp_out = (ocp_qp_out *) ssGetPWork(S)[2];
    ocp_qp_solver *qp_solver = (ocp_qp_solver *) ssGetPWork(S)[4];

    const double *x0 = ssGetInputPortRealSignal(S, 0);
    d_cvt_colmaj_to_ocp_qp_lbx(0, x0, qp_in);
    d_cvt_colmaj_to_ocp_qp_ubx(0, x0, qp_in);

    int status = ocp_qp_solve(qp_solver, qp_in, qp_out);

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;

    double *u0_opt = ssGetOutputPortRealSignal(S, 0);
    double *time = ssGetOutputPortRealSignal(S, 1);
    double *status_out = ssGetOutputPortRealSignal(S, 2);

    d_cvt_ocp_qp_sol_to_colmaj_u(qp_out, u0_opt, 0);
    *time = info->total_time;
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

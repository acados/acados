#define S_FUNCTION_NAME   acados_solver_sfunction_{{ra.model_name}}
#define S_FUNCTION_LEVEL  2

#define MDL_START

// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// TODO(oj): remove, when setters for Cyt,idxb available
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// example specific
#include "{{ ra.model_name }}_model/{{ ra.model_name }}_model.h"
#include "acados_solver_{{ ra.model_name }}.h"

#include "simstruc.h"

#define SAMPLINGTIME -1

static void mdlInitializeSizes (SimStruct *S)
{
    // specify the number of continuous and discrete states 
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    // specify the number of input ports 
    if ( !ssSetNumInputPorts(S, 1) )
        return;

    // specify the number of output ports 
    if ( !ssSetNumOutputPorts(S, 4) )
        return;

    // specify dimension information for the input ports 
    ssSetInputPortVectorDimension(S, 0, {{ ra.dims.nx }});

    // specify dimension information for the output ports 
    ssSetOutputPortVectorDimension(S, 0, {{ ra.dims.nu }} ); // optimal input
    ssSetOutputPortVectorDimension(S, 1, 1 );                // solver status
    ssSetOutputPortVectorDimension(S, 2, 1 );                // KKT residuals
    ssSetOutputPortVectorDimension(S, 3, {{ ra.dims.nx }} ); // first state

    // specify the direct feedthrough status 
    ssSetInputPortDirectFeedThrough(S, 0, 1); // current state x0

    // one sample time 
    ssSetNumSampleTimes(S, 1);
    }


#if defined(MATLAB_MEX_FILE)

#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetInputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetOutputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

    #endif /* MATLAB_MEX_FILE */


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLINGTIME);
    ssSetOffsetTime(S, 0, 0.0);
}


static void mdlStart(SimStruct *S)
{
    acados_create();
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    // get pointers to acados internal structures
    ocp_nlp_in * _nlp_in =  acados_get_nlp_in();
    ocp_nlp_out * _nlp_out =  acados_get_nlp_out();
    ocp_nlp_solver * _nlp_solver =  acados_get_nlp_solver();
    void * _nlp_opts =  acados_get_nlp_opts();
    ocp_nlp_solver_config * _nlp_config =  acados_get_nlp_config();
    ocp_nlp_dims * _nlp_dims =  acados_get_nlp_dims();

    // get input signals
    InputRealPtrsType in_x0_sign;
    
    // local buffers
    real_t in_x0[{{ ra.dims.nx }}];

    in_x0_sign = ssGetInputPortRealSignalPtrs(S, 0);

    // copy signals into local buffers
    for (int i = 0; i < {{ ra.dims.nx }}; i++) in_x0[i] = (double)(*in_x0_sign[i]);

    for (int i = 0; i < 4; i++) ssPrintf("x0[%d] = %f\n", i, in_x0[i]);
    ssPrintf("\n");

    // set initial condition
    ocp_nlp_constraints_bounds_set(_nlp_config, _nlp_dims, _nlp_in, 0, "lbx", in_x0);
    ocp_nlp_constraints_bounds_set(_nlp_config, _nlp_dims, _nlp_in, 0, "ubx", in_x0);
    
    // assign pointers to output signals 
    real_t *out_u0, *out_status, *out_KKT_res, *out_x1;

    out_u0      = ssGetOutputPortRealSignal(S, 0);
    out_status  = ssGetOutputPortRealSignal(S, 1);
    out_KKT_res = ssGetOutputPortRealSignal(S, 2);
    out_x1      = ssGetOutputPortRealSignal(S, 3);
    
    // get pointers to acados structures 
    _nlp_opts = acados_get_nlp_opts();
    _nlp_dims = acados_get_nlp_dims();
    _nlp_config = acados_get_nlp_config();
    _nlp_out = acados_get_nlp_out();
    _nlp_in = acados_get_nlp_in();

    // call acados_solve()
    acados_solve();
    
    // get solution
    ssPrintf("%p\n",(void*)_nlp_config);
    ssPrintf("%p\n",(void*)_nlp_dims);
    ssPrintf("%p\n",(void*)_nlp_out);
    ocp_nlp_out_get(_nlp_config, _nlp_dims, _nlp_out, 0, "u", (void *) out_u0);

    // get next state
    ocp_nlp_out_get(_nlp_config, _nlp_dims, _nlp_out, 1, "x", (void *) out_x1);

}

static void mdlTerminate(SimStruct *S)
{
    acados_free();
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif

#define S_FUNCTION_NAME   acado_solver_sfun
#define S_FUNCTION_LEVEL  2

#define MDL_START

#define VARYING_REFERENCE 1
#define VARYING_DATA 1
#define RESET_WHEN_NAN 0   // no reset ( = 0 ), reset with previous traj ( = 1 ), reset to initial values ( = 2 )
#define KKT_MIN -1
#define KKT_MAX 1e15

#include "acado_common.h"
#include "acado_auxiliary_functions.h"
#include "simstruc.h"

#define SAMPLINGTIME -1

static void mdlInitializeSizes (SimStruct *S)
{
    // specify the number of continuous and discrete states 
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    // specify the number of input ports 
    // state of the system: x_0 
    if ( !ssSetNumInputPorts(S, 1) )
        return;

    // specify the number of output ports 
    // optimal input: u_0^*
    // solver status: status
    // KKT residuals: KKT_res
    
    if ( !ssSetNumOutputPorts(S, 3) )
        return;

    // // specify the number of parameters 
    // ssSetNumSFcnParams(S, 5);
    // if ( ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S) )
    //     return;

    // specify dimension information for the input ports 
    ssSetInputPortVectorDimension(S, 0, {{ ra.dims.nx }});
    // #if VARYING_REFERENCE
    // ssSetInputPortVectorDimension(S, 1, ACADO_N*ACADO_NY+ACADO_NYN);
    // #else
    ssSetInputPortVectorDimension(S, 1, ACADO_NY);
    // #endif
    ssSetInputPortVectorDimension(S, 2, ACADO_NY*ACADO_NY);
    ssSetInputPortVectorDimension(S, 3, ACADO_NYN*ACADO_NYN);
    ssSetInputPortVectorDimension(S, 4, 1);
    ssSetInputPortVectorDimension(S, 5, 1);
    ssSetInputPortVectorDimension(S, 6, (ACADO_N+1)*ACADO_NOD);

    /* Specify dimension information for the output ports */
    ssSetOutputPortVectorDimension(S, 0, ACADO_NU );                // first input
    ssSetOutputPortMatrixDimensions(S, 1, ACADO_N+1, ACADO_NX );    // predicted state trajectory
    ssSetOutputPortMatrixDimensions(S, 2, ACADO_N, ACADO_NU );      // optimal input trajectory
    ssSetOutputPortVectorDimension(S, 3, 1 );                       //
    ssSetOutputPortVectorDimension(S, 4, 1 );
    ssSetOutputPortVectorDimension(S, 5, 1 );
    ssSetOutputPortVectorDimension(S, 6, 1 );
    ssSetOutputPortVectorDimension(S, 7, 1 );

    /* Specify the direct feedthrough status */
    ssSetInputPortDirectFeedThrough(S, 0, 1); // current state x0
    ssSetInputPortDirectFeedThrough(S, 1, 1); // reference
    ssSetInputPortDirectFeedThrough(S, 2, 1); // weghting matrix
    ssSetInputPortDirectFeedThrough(S, 3, 1); // terminal weghting matrix
    ssSetInputPortDirectFeedThrough(S, 4, 1); // number of sqp iterations
    ssSetInputPortDirectFeedThrough(S, 5, 1); // shifting
    ssSetInputPortDirectFeedThrough(S, 6, 1); // online data (m_load)


    /* One sample time */
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
    int i, j, k;

    InputRealPtrsType in_ref, in_od;
    double *xInit, *zInit, *uInit, *Smat, *SNmat;

    /* get inputs and perform feedback step */
    in_ref = ssGetInputPortRealSignalPtrs(S, 1);
    in_od = ssGetInputPortRealSignalPtrs(S, 6);


    xInit = mxGetPr( ssGetSFcnParam(S, 0) );
    zInit = mxGetPr( ssGetSFcnParam(S, 1) );
    uInit = mxGetPr( ssGetSFcnParam(S, 2) );


    for( i=0; i < ACADO_N+1; ++i ) {
        for( j=0; j < ACADO_NX; ++j )acadoVariables.x[i*ACADO_NX+j] = xInit[j];
    }

    for( i=0; i < ACADO_N; ++i ) {
        for( j=0; j < ACADO_NXA; ++j ) acadoVariables.z[i*ACADO_NXA+j] = zInit[j];
    }

    for( i=0; i < ACADO_N; ++i ) {
        for( j=0; j < ACADO_NU; ++j ) acadoVariables.u[i*ACADO_NU+j] = uInit[j];
    }

    for( i=0; i < ACADO_N; ++i ) {
        for( j=0; j < ACADO_NY; ++j ) acadoVariables.y[i*ACADO_NY+j] = (double)(*in_ref[i*ACADO_NY+j]);
    }
    for( i=0; i < ACADO_NYN; ++i ) acadoVariables.yN[i] = (double)(*in_ref[ACADO_N*ACADO_NY+i]);

    Smat = mxGetPr( ssGetSFcnParam(S, 3) );
    SNmat = mxGetPr( ssGetSFcnParam(S, 4) );


    for( i = 0; i < (ACADO_NY); ++i )  {
        for( j = 0; j < ACADO_NY; ++j ) {
            acadoVariables.W[i*ACADO_NY+j] = Smat[i*ACADO_NY+j];
        }
    }

    for( i = 0; i < (ACADO_NYN); ++i )  {
        for( j = 0; j < ACADO_NYN; ++j ) {
            acadoVariables.WN[i*ACADO_NYN+j] = SNmat[i*ACADO_NYN+j];
        }
    }

    for( i=0; i < ACADO_N+1; ++i ) {
        for( j=0; j < ACADO_NOD; ++j ) acadoVariables.od[i*ACADO_NOD+j] = (double)(*in_od[i*ACADO_NOD+j]);
    }

    acado_initializeSolver();

}

extern int acado_getNWSR();

static void mdlOutputs(SimStruct *S, int_T tid)
{
    int i, j, status, sumIter;
    real_t kkt;
	real_t* xEnd = NULL;
	real_t* uEnd = NULL;

    InputRealPtrsType in_x, in_ref, in_W, in_WN, in_bValues, in_nIter, in_shifting, in_od;
    real_t *out_u0, *out_xTraj, *out_uTraj, *out_kktTol, *out_cpuTime, *out_status, *out_nIter, *out_objVal, *xInit, *zInit, *uInit;

    out_u0     = ssGetOutputPortRealSignal(S, 0);
    out_xTraj  = ssGetOutputPortRealSignal(S, 1);
    out_uTraj  = ssGetOutputPortRealSignal(S, 2);
    out_status = ssGetOutputPortRealSignal(S, 3);
    out_kktTol = ssGetOutputPortRealSignal(S, 4);
    out_cpuTime = ssGetOutputPortRealSignal(S, 5);
    out_nIter = ssGetOutputPortRealSignal(S, 6);
    out_objVal = ssGetOutputPortRealSignal(S, 7);

    /* get inputs and perform feedback step */
    in_x    = ssGetInputPortRealSignalPtrs(S, 0);
    in_ref  = ssGetInputPortRealSignalPtrs(S, 1);
    in_W  = ssGetInputPortRealSignalPtrs(S, 2);
    in_WN  = ssGetInputPortRealSignalPtrs(S, 3);
    in_nIter   = ssGetInputPortRealSignalPtrs(S, 4);
    in_shifting   = ssGetInputPortRealSignalPtrs(S, 5);
    in_od   = ssGetInputPortRealSignalPtrs(S, 6);



    for( i=0; i < ACADO_NX; ++i ) {
        acadoVariables.x0[i] = (double)(*in_x[i]);
        // printf("x0[%i] = %f\n",i,(double)(*in_x[i]));
    };

    for( i=0; i < ACADO_N; ++i ) {
        for( j=0; j < ACADO_NY; ++j ) {
            acadoVariables.y[i*ACADO_NY+j] = (double)(*in_ref[i*ACADO_NY+j]);
            // printf("ref y[%i] = %f\n",i*ACADO_NY+j,(double)(*in_ref[i*ACADO_NY+j]));
    }
    }

    for( i=0; i < ACADO_NYN; ++i ){
        acadoVariables.yN[i] = (double)(*in_ref[ACADO_N*ACADO_NY+i]);
        // printf("ref yN[%i] = %f\n",ACADO_N*ACADO_NY+i,(double)(*in_ref[ACADO_N*ACADO_NY+i]));
    }
//
    for( i=0; i < ACADO_NY; ++i ) {
        for( j=0; j < ACADO_NY; ++j ){
            acadoVariables.W[i*ACADO_NY+j] = (double)(*in_W[i*ACADO_NY+j]);
            // printf("ref W[%i] = %f\n",i*ACADO_NY+j,(double)(*in_W[i*ACADO_NY+j]));

        }
    }

    for( i=0; i < ACADO_NYN; ++i ) {
        for( j=0; j < ACADO_NYN; ++j ) {
            acadoVariables.WN[i*ACADO_NYN+j] = (double)(*in_WN[i*ACADO_NYN+j]);
            // printf("ref WN[%i] = %f\n",i*ACADO_NY+j,(double)(*in_W[i*ACADO_NY+j]));

        }
    }

    for( i=0; i < ACADO_N+1; ++i ) {
        for( j=0; j < ACADO_NOD; ++j ) acadoVariables.od[i*ACADO_NOD+j] = (double)(*in_od[i*ACADO_NOD+j]);
    }

    // printf("before acado_preparationStep");
    acado_preparationStep( );

    status = acado_feedbackStep(  );
    // printf("before acado_feedbackStep");

    sumIter = (int)acado_getNWSR();
    // printf("before acado_getNWSR");

//     if(in_shifting)
//     {
//         acado_shiftStates(1, 0, 0);
//         acado_shiftControls(0);
//     }
//
//     for( i=0; i < ACADO_N+1; ++i ) {
//         for( j = 0; j < ACADO_NX; ++j ) {
//             out_xTraj[j*(ACADO_N+1)+i] = acadoVariables.x[i*ACADO_NX+j];
//         }
//     }
//     for( i=0; i < ACADO_N; ++i ) {
//         for( j = 0; j < ACADO_NU; ++j ) {
//             out_uTraj[j*(ACADO_N)+i] = acadoVariables.u[i*ACADO_NU+j];
//         }
//     }

    /* return outputs*/
    for( i=0; i < ACADO_NU; ++i ) out_u0[i] = acadoVariables.u[i];

    for( i=0; i < ACADO_N+1; ++i ) {
        for( j = 0; j < ACADO_NX; ++j ) {
            out_xTraj[j*(ACADO_N+1)+i] = acadoVariables.x[i*ACADO_NX+j];
        }
    }
    for( i=0; i < ACADO_N; ++i ) {
        for( j = 0; j < ACADO_NU; ++j ) {
            out_uTraj[j*(ACADO_N)+i] = acadoVariables.u[i*ACADO_NU+j];
        }
    }

    // printf("before acado_getKKT");
    out_kktTol[0] = acado_getKKT( );
    // printf("done with acado_getKKT");
    out_status[0] = status;
    out_cpuTime[0] = -1;
    out_nIter[0] = (double) sumIter;
    // printf("before acado_getObjective");
    out_objVal[0] = acado_getObjective( );
    // printf("after acado_getObjective");

}


// #define MDL_UPDATE
// static void mdlUpdate(SimStruct *S, int_T tid)
// {
// //     timer t;
//
// //     tic( &t );
//     acado_preparationStep( );
// //     timePrep = toc( &t );
//
// //     printf("Timing RTI iteration (MPC preparation):                   %.3g ms   \n", timePrep*1e3);
// //     printf("---------------------------------------------------------------\n");
// }


static void mdlTerminate(SimStruct *S)
{
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif

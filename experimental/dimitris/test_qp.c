#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "acados/utils/types.h"
#include "acados/utils/timing.h"
#include "acados/ocp_qp/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "test/test_utils/read_ocp_qp_in.h"

#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#define OOQP_WORK 2  // 1: obsolete, 2: chunk of memory

int_t main( ) {
    /* code */
    int N, return_value;

    #if SOLVER == 1
    ocp_qp_ooqp_args args;
    ocp_qp_ooqp_memory mem;
    #if OOQP_WORK == 1
    ocp_qp_ooqp_workspace work;
    #elif OOQP_WORK == 2
    char *work;
    #endif
    args.printLevel = 0;
    args.fixHessian = 0;
    args.fixHessianSparsity = 0;
    args.fixDynamics = 0;
    args.fixDynamicsSparsity = 0;
    args.fixInequalities = 0;
    args.fixInequalitiesSparsity = 0;
    #elif SOLVER == 2
    ocp_qp_condensing_qpoases_args args;
    args.dummy = 42;
    #elif SOLVER == 3
    ocp_qp_hpmpc_opts args;
    args.tol = 1e-8;
    args.max_iter = 20;
    args.mu0 = 0.0;
    args.warm_start = 0;
    args.N2 = N;
    double inf_norm_res[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    args.inf_norm_res = &inf_norm_res[0];
    #endif

    int_t BOUNDS = 1;
    int_t CONSTRAINTS = 1;
    int_t MPC = 1;
    int_t QUIET = 1;

    ocp_qp_in *qp_in = read_ocp_qp_in("LTV/", BOUNDS, CONSTRAINTS, MPC, QUIET);
    ocp_qp_out *qp_out = ocp_qp_out_create(qp_in->N, (int*)qp_in->nx, (int*)qp_in->nu, (int*)qp_in->nb, (int*)qp_in->nc);

    N = qp_in->N;
    int nPrimalVars = 0;
    for (int kk = 0; kk < N; kk++) {
        nPrimalVars += qp_in->nx[kk] + qp_in->nu[kk];
    }
    nPrimalVars += qp_in->nx[N];

    #if SOLVER == 1
    ocp_qp_ooqp_create_memory(qp_in, &args, &mem);
    #if OOQP_WORK == 1
    ocp_qp_ooqp_create_workspace(qp_in, &args, &work);
    #elif OOQP_WORK == 2
    int_t work_space_size = ocp_qp_ooqp_calculate_workspace_size(qp_in, &args);
    printf("\nwork space size: %d bytes\n", work_space_size);
    work = (void*)malloc(work_space_size);
    #endif
    #elif SOLVER == 2
    initialise_qpoases(qp_in);
    #elif SOLVER == 3
    int work_space_size =
        ocp_qp_hpmpc_workspace_size_bytes(N, (int_t *)qp_in->nx, (int_t *)qp_in->nu, (int_t *)qp_in->nb, (int_t *)qp_in->nc, (int_t **)qp_in->idxb, &args);
    printf("\nwork space size: %d bytes\n", work_space_size);
    void *work = (void*)malloc(work_space_size);
    #endif

    #if SOLVER == 1
    #if OOQP_WORK == 1
    return_value = ocp_qp_ooqp(qp_in, qp_out, &args, &mem, &work);
    #elif OOQP_WORK == 2
    return_value = ocp_qp_ooqp(qp_in, qp_out, &args, &mem, work);
    return_value = ocp_qp_ooqp(qp_in, qp_out, &args, &mem, work);
    #endif
    #elif SOLVER == 2
    return_value = ocp_qp_condensing_qpoases(qp_in, qp_out, &args, NULL);
    #elif SOLVER ==3
    return_value = ocp_qp_hpmpc(qp_in, qp_out, &args, work);
    #endif

    printf("\nRETURN VALUE = %d LAST ELEMENT OF SOLUTION = %f\n\n", return_value,qp_out->x[0][nPrimalVars-1]);
    if (!QUIET) {
        printf("\nSOLUTION= \n");
        for (int ii=0; ii <nPrimalVars; ii++) printf("%10.8f\n",qp_out->x[0][ii]);
    }

    #if SOLVER == 1
    ocp_qp_ooqp_free_memory(&mem);
    #if OOQP_WORK == 1
    ocp_qp_ooqp_free_workspace(&work);
    #elif OOQP_WORK == 2
    free(work);
    #endif
    #elif SOLVER == 3
    free(work);
    #endif

    free(qp_in);
    free(qp_out);

    return 0;
}

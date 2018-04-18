#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>

#include "acados/utils/types.h"
#include "acados/utils/timing.h"
#include "acados/utils/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "test/test_utils/read_ocp_qp_in.h"

#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"

#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#define TOL 1e-6

real_t error_in_primal_solution(int_t n, real_t *v1, real_t *v2) {
    real_t error = 0;
    int_t i;

    // calculate 2-norm
    for (i = 0; i < n; i++) {
        error += pow(v1[i] - v2[i], 2);
    }
    error = sqrt(error);
    printf(">>>>>>>>>>> ERROR = %e\n\n", error);
    return error;
}

int_t main( ) {
    /* code */
    int_t N, return_value;
    int_t iScenario, iSolver, iProblem, nScenarios, nSolvers, nProblems;
    int_t BOUNDS, CONSTRAINTS, MPC, QUIET;
    int_t work_space_size, ooqp_work_space_size, hpmpc_work_space_size;
    real_t *sol;
    char fname[256];

    // define input-output
    ocp_qp_in qp_in;
    ocp_qp_out qp_out;

    // define arguments for all solvers
    ocp_qp_ooqp_args ooqp_args;
    ocp_qp_condensing_qpoases_args qpoases_args;
    ocp_qp_hpmpc_opts hpmpc_args;

    ooqp_args.printLevel = 0;
    ooqp_args.fixHessian = 0;
    ooqp_args.fixHessianSparsity = 0;
    ooqp_args.fixDynamics = 0;
    ooqp_args.fixDynamicsSparsity = 0;
    ooqp_args.fixInequalities = 0;
    ooqp_args.fixInequalitiesSparsity = 0;

    qpoases_args.dummy = 42;

    hpmpc_args.tol = 1e-15;
    hpmpc_args.max_iter = 20;
    hpmpc_args.mu0 = 0.0;
    hpmpc_args.warm_start = 0;
    double inf_norm_res[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    hpmpc_args.inf_norm_res = &inf_norm_res[0];

    // define memory for all solvers (that have it implemented..)
    ocp_qp_ooqp_memory ooqp_mem;

    void *work;

    const char *problems[] = {"LTI", "LTV"};
    const char *scenarios[] = {"UNCONSTRAINED", "ONLY_BOUNDS", "ONLY_AFFINE", "CONSTRAINED"};
    const char *solvers[] = {"qpoases", "ooqp", "hpmpc"};

    nScenarios = (int_t)(sizeof(scenarios)/sizeof(scenarios[0]));
    nProblems = (int_t)(sizeof(problems)/sizeof(problems[0]));
    nSolvers = (int_t)(sizeof(solvers)/sizeof(solvers[0]));

    QUIET = 1;
    MPC = 1;

    for (iScenario = 0; iScenario < nScenarios; iScenario++) {
        for (iProblem = 0; iProblem < nProblems; iProblem++) {
            if (strcmp(scenarios[iScenario], "UNCONSTRAINED") == 0) {
                BOUNDS = 0; CONSTRAINTS = 0;
                snprintf(fname, sizeof(fname), "%s%s", problems[iProblem], "/sol_only_x0.txt");
            } else if (strcmp(scenarios[iScenario], "ONLY_BOUNDS") == 0) {
                BOUNDS = 1; CONSTRAINTS = 0;
                snprintf(fname, sizeof(fname), "%s%s", problems[iProblem], "/sol_only_bounds.txt");
            } else if (strcmp(scenarios[iScenario], "ONLY_AFFINE") == 0) {
                BOUNDS = 0; CONSTRAINTS = 1;
                snprintf(fname, sizeof(fname), "%s%s", problems[iProblem], "/sol_only_ineq.txt");
            } else {
                BOUNDS = 1; CONSTRAINTS = 1;
                snprintf(fname, sizeof(fname), "%s%s", problems[iProblem], "/sol_constrained.txt");
            }

            int_t NEW_PROBLEM = 1;

            for (iSolver = 0; iSolver < nSolvers; iSolver++) {
                printf(">>>>> SOLVING %s QP (%s) with %s\n", problems[iProblem], scenarios[iScenario], solvers[iSolver]);

                read_ocp_qp_in(&qp_in,  problems[iProblem], BOUNDS, CONSTRAINTS, MPC, QUIET);
                allocate_ocp_qp_out(&qp_in, &qp_out);

                // calculate number of primal variables
                N = qp_in.N;
                int nPrimalVars = 0;
                for (int kk = 0; kk < N; kk++) {
                    nPrimalVars += qp_in.nx[kk] + qp_in.nu[kk];
                }
                nPrimalVars += qp_in.nx[N];
                sol = (real_t*)malloc(sizeof(real_t)*nPrimalVars);

                hpmpc_args.N2 = N;

                // printf(">>>>>>>>>>> READING SOLUTION FROM %s\n", fname);
                read_double_vector_from_txt(sol, nPrimalVars, fname);

                // calculate workspace size for hpmpc and ooqp for the current problem (once) and allocate maximum
                if (NEW_PROBLEM) {
                    printf("\n ----- ALLOCATING COMMON WORKSPACE FOR CURRENT PROBLEM ----- \n");
                    // printf("calculating OOQP workspace size... \n");
                    ooqp_work_space_size = ocp_qp_ooqp_calculate_workspace_size(&qp_in, &ooqp_args);
                    // printf("calculating HPMPC workspace size... \n");
                    hpmpc_work_space_size = ocp_qp_hpmpc_workspace_size_bytes(N, (int_t *)qp_in.nx,
                    (int_t *)qp_in.nu, (int_t *)qp_in.nb, (int_t *)qp_in.nc, (int_t **)qp_in.idxb,
                    &hpmpc_args);
                    work_space_size = max(ooqp_work_space_size, hpmpc_work_space_size);
                    // printf("ALLOCATING %d bytes max(%d, %d)\n", work_space_size, ooqp_work_space_size, hpmpc_work_space_size);
                    work = (void*)malloc(work_space_size);
                    NEW_PROBLEM = 0;
                }

                if (strcmp(solvers[iSolver], "qpoases") == 0) {
                    initialise_qpoases(&qp_in);
                    return_value = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &qpoases_args, NULL);

                } else if (strcmp(solvers[iSolver], "ooqp") == 0) {
                    ocp_qp_ooqp_create_memory(&qp_in, &ooqp_args, &ooqp_mem);
                    return_value = ocp_qp_ooqp(&qp_in, &qp_out, &ooqp_args, &ooqp_mem, work);
                    ocp_qp_ooqp_free_memory(&ooqp_mem);
                } else if (strcmp(solvers[iSolver], "hpmpc") == 0) {
                    return_value = ocp_qp_hpmpc(&qp_in, &qp_out, &hpmpc_args, work);
                }

                printf("\n>>>>>>>>>>> RETURN VALUE = %d\n\n", return_value);
                if (!QUIET) {
                    printf("\nSOLUTION= \n");
                    for (int ii=0; ii <nPrimalVars; ii++) printf("%10.8f\n",qp_out.x[0][ii]);
                }

                assert(return_value == 0);
                assert(error_in_primal_solution(nPrimalVars, qp_out.x[0], sol) < TOL);

                free_ocp_qp_in(&qp_in);
                free_ocp_qp_out(&qp_out);
                free(sol);
            }
            free(work);
        }
    }
    printf(" >>>>>>>>>>>>>> AWESOME! <<<<<<<<<<<<<<<<\n\n");
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <qpDUNES.h>

#include "acados/utils/types.h"
#include "acados/utils/timing.h"
#include "acados/utils/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "test/test_utils/read_ocp_qp_in.h"

#define INFTY 1.0e12

int main(int argc, char const *argv[]) {
    /* general vars */
    int_t iter, ii, kk, nx, nu, N;
    uint_t *nD;

    /* qpDUNES vars */
    return_t statusFlag;
    qpData_t qpData;
    boolean_t isLTI;
    qpOptions_t qpOptions = qpDUNES_setupDefaultOptions();
    real_t *zLow, *zUpp, *g;

    ocp_qp_qpdunes_args qpdunes_args;
    ocp_qp_qpdunes_memory qpdunes_mem;

    /* acados vars */
    ocp_qp_in qp_in;
    ocp_qp_out qp_out;

    int_t BOUNDS = 1;
    int_t CONSTRAINTS = 0;  // always for qpDUNES with clipping
    int_t MPC = 1;
    int_t QUIET = 1;

    char *test_problem = "LTI_q0";

    /* read and allocate data */
    read_ocp_qp_in(&qp_in, test_problem, BOUNDS, CONSTRAINTS, MPC, QUIET);
    allocate_ocp_qp_out(&qp_in, &qp_out);

    N  = qp_in.N;
    nx = qp_in.nx[0];
    nu = qp_in.nu[0];

    qpdunes_args.options = qpDUNES_setupDefaultOptions();
    ocp_qp_qpdunes_create_memory(&qp_in, &qpdunes_args, &qpdunes_mem);

    qpData = qpdunes_mem.qpData;

    /* solve problem */
    printf("BEFORE SOLVE\n---------------------------------------------------------");
	statusFlag = qpDUNES_solve(&(qpdunes_mem.qpData));
	if (statusFlag != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND)
	{
		printf("qpDUNES failed to solve the QP. Error code: %d\n", statusFlag);
		return (int)statusFlag;
	}

    /* write out solution */
    // qpDUNES_getPrimalSol(&qpData, zOpt);

    printf("\nz_opt (qpDUNES) = \n");
    for(kk = 0; kk < N; kk++) {
        for (ii = 0; ii < nx+nu; ii++) printf("%5.3f\n", qpData.intervals[kk]->z.data[ii]);
    }
    for (ii = 0; ii < nx; ii++) printf("%5.3f\n", qpData.intervals[kk]->z.data[ii]);
    printf("\n");

	qpDUNES_cleanup( &qpData );

    // qpDUNES_printf("Default iterations: %d", qpOptions.maxIter);
    free(zLow); free(zUpp); free(g);

    // --------------> SOLVE WITH OOQP TO COMPARE
    read_ocp_qp_in(&qp_in, test_problem, BOUNDS, CONSTRAINTS, MPC, QUIET);

    ocp_qp_ooqp_args ooqp_args;
    ocp_qp_ooqp_memory ooqp_mem;

    ooqp_args.workspaceMode = 2;
    ooqp_args.printLevel = 0;

    ocp_qp_ooqp_create_memory(&qp_in, &ooqp_args, &ooqp_mem);
    int_t work_space_size = ocp_qp_ooqp_calculate_workspace_size(&qp_in, &ooqp_args);
    printf("\nwork space size: %d bytes\n", work_space_size);
    void *work = (void*)malloc(work_space_size);
    statusFlag = ocp_qp_ooqp(&qp_in, &qp_out, &ooqp_args, &ooqp_mem, work);
    ocp_qp_ooqp_free_memory(&ooqp_mem);
    free(work);

    printf("\nz_opt (OOQP) (FLAG = %d) = \n", statusFlag);
    for (ii = 0; ii < N*(nx+nu)+nx; ii++) printf("%5.3f\n",qp_out.x[0][ii]);

    free_ocp_qp_in(&qp_in);
    free_ocp_qp_out(&qp_out);
    return 0;
}

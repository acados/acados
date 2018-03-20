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
    boolean_t isLTI;

    /* acados vars */
    ocp_qp_in qp_in;
    ocp_qp_out qp_out_qpdunes, qp_out_ooqp;

    ocp_qp_qpdunes_opts qpdunes_args;
    ocp_qp_qpdunes_memory qpdunes_mem;
    int_t qpdunes_workspace_size;
    void *qpdunes_work;

    /* sim vars */
    int_t BOUNDS = 1;
    int_t CONSTRAINTS = 1;
    int_t MPC = 1;
    int_t QUIET = 1;

    char *test_problem = "LTV/";

    /* read and allocate data */
    read_ocp_qp_in(&qp_in, test_problem, BOUNDS, CONSTRAINTS, MPC, QUIET);
    allocate_ocp_qp_out(&qp_in, &qp_out_qpdunes);
    allocate_ocp_qp_out(&qp_in, &qp_out_ooqp);

    N  = qp_in.N;
    nx = qp_in.nx[0];
    nu = qp_in.nu[0];

    /* setup QP */
    ocp_qp_qpdunes_create_arguments(&qpdunes_args, QPDUNES_DEFAULT_ARGUMENTS);

    qpdunes_workspace_size = ocp_qp_qpdunes_calculate_workspace_size(&qp_in, &qpdunes_args);
    qpdunes_work = (void*)malloc(qpdunes_workspace_size);

    ocp_qp_qpdunes_create_memory(&qp_in, &qpdunes_args, &qpdunes_mem);

    ocp_qp_qpdunes(&qp_in, &qp_out_qpdunes, &qpdunes_args, &qpdunes_mem, qpdunes_work);

    ocp_qp_qpdunes_free_memory(&qpdunes_mem);
    free(qpdunes_work);

    printf("\nz_opt (qpDUNES) = \n");
    for (ii = 0; ii < N*(nx+nu)+nx; ii++) printf("%5.3f\n",qp_out_qpdunes.x[0][ii]);

    // --------------> SOLVE WITH OOQP TO COMPARE

    ocp_qp_ooqp_args ooqp_args;
    ocp_qp_ooqp_memory ooqp_mem;

    ooqp_args.printLevel = 0;

    ocp_qp_ooqp_create_memory(&qp_in, &ooqp_args, &ooqp_mem);
    int_t work_space_size = ocp_qp_ooqp_calculate_workspace_size(&qp_in, &ooqp_args);
    printf("\nwork space size: %d bytes\n", work_space_size);
    void *work = (void*)malloc(work_space_size);
    statusFlag = ocp_qp_ooqp(&qp_in, &qp_out_ooqp, &ooqp_args, &ooqp_mem, work);
    ocp_qp_ooqp_free_memory(&ooqp_mem);
    free(work);

    printf("\nz_opt (OOQP) (FLAG = %d) = \n", statusFlag);
    for (ii = 0; ii < N*(nx+nu)+nx; ii++) printf("%5.3f\n",qp_out_ooqp.x[0][ii]);

    printf("\nERROR = \n");
    for (ii = 0; ii < N*(nx+nu)+nx; ii++) printf("%10.9f\n",qp_out_qpdunes.x[0][ii]-qp_out_ooqp.x[0][ii]);


    free_ocp_qp_in(&qp_in);
    free_ocp_qp_out(&qp_out_qpdunes);
    free_ocp_qp_out(&qp_out_ooqp);
    return 0;
}

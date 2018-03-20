#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/print.h"
#include "acados_c/dense_qp_interface.h"

#define QP_HORIZON 5

int main() {

    double H[] = {3., 1., 1., 4.};
    double g[] = {2., 1.};

    int idxb[] = {0};
    double d_lb[] = {-1.};
    double d_ub[] = {1.};

    double C[] = {1.2, 1.3};
    double d_lg[] = {-2.};
    double d_ug[] = {2.};

    dense_qp_solver_plan plan;
    plan.qp_solver = DENSE_QP_HPIPM;

    qp_solver_config *config = dense_qp_config_create(&plan);

    dense_qp_dims dims;
    dims.nv = 2;
    dims.ne = 0;  // TODO(dimitris): is this even supported in acados?
    dims.nb = 1;
    dims.ng = 1;
    dims.ns = 0;

    dense_qp_in *qp_in = dense_qp_in_create(config, &dims);

    d_cvt_colmaj_to_dense_qp(H, g, NULL, NULL, idxb, d_lb, d_ub, C, d_lg, d_ug, NULL, NULL, NULL, NULL, NULL, qp_in);

    print_dense_qp_in(qp_in);

    void *opts = dense_qp_opts_create(config, &dims);

    dense_qp_out *qp_out = dense_qp_out_create(config, &dims);

    dense_qp_solver *qp_solver = dense_qp_create(config, &dims, opts);

    int acados_return = dense_qp_solve(qp_solver, qp_in, qp_out);

    if (acados_return != ACADOS_SUCCESS)
        return -1;

    // TODO(dimitris): implement this
    // print_dense_qp_out(qp_out);

    /************************************************
     * compute inf norm of residuals
     ************************************************/

    double res[4];
    dense_qp_inf_norm_residuals(&dims, qp_in, qp_out, res);
    printf("\ninf norm res: %e, %e, %e, %e\n\n", res[0], res[1], res[2], res[3]);

    free(config);
    free(opts);
    free(qp_in);
    free(qp_out);
    free(qp_solver);
}

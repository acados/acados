
#include <stdio.h>

#include <stdlib.h>
#include <stddef.h>

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/utils/print.h"
#include "acados/utils/mem.h"
#include "acados_c/dense_qp_interface.h"
// hpipm
#include "hpipm/include/hpipm_d_dense_qp_kkt.h"

static dense_qp_res *dense_qp_res_create(dense_qp_dims *dims)
{
    int size = dense_qp_res_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_res *qp_res = dense_qp_res_assign(dims, ptr);
    return qp_res;
}

static dense_qp_res_ws *dense_qp_res_workspace_create(dense_qp_dims *dims)
{
    int size = dense_qp_res_workspace_calculate_size(dims);
    void *ptr = acados_malloc(size, 1);
    dense_qp_res_ws *res_ws = dense_qp_res_workspace_assign(dims, ptr);
    return res_ws;
}

int main() {

    double H[] = {3., 1., 1., 4.};
    double g[] = {2., 1.};

    int idxb[] = {0};
    double d_lb[] = {-1.};
    double d_ub[] = {1.};

    double C[] = {1.2, 1.3};
    double d_lg[] = {-2.};
    double d_ug[] = {2.};

    int idxs[] = {0, 1};
    double d_ls[] = {0.0, 0.0};
    double d_us[] = {0.0, 0.0};

    double Zl[] = {1.0, 1.0};
    double zl[] = {1.0, 1.0};
    double Zu[] = {1.0, 1.0};
    double zu[] = {1.0, 1.0};

    dense_qp_solver_plan plan;
    plan.qp_solver = DENSE_QP_HPIPM;

    qp_solver_config *config = dense_qp_config_create(&plan);

    dense_qp_dims dims;
    dims.nv = 2;
    dims.ne = 0;  // TODO(dimitris): is this even supported in acados?
    dims.nb = 1;
    dims.ng = 1;
    dims.ns = 2;
//    dims.nsb = 1;
//    dims.nsg = 1;

    dense_qp_in *qp_in = dense_qp_in_create(config, &dims);

    d_cvt_colmaj_to_dense_qp(H, g, NULL, NULL, idxb, d_lb, d_ub, C, d_lg, d_ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, qp_in);

    print_dense_qp_in(qp_in);

    void *opts = dense_qp_opts_create(config, &dims);

    dense_qp_out *qp_out = dense_qp_out_create(config, &dims);

    dense_qp_solver *qp_solver = dense_qp_create(config, &dims, opts);

	// overwrite default opts
	dense_qp_hpipm_opts *hpipm_opts = opts;
	hpipm_opts->hpipm_opts->mu0 = 1e1;

    int acados_return = dense_qp_solve(qp_solver, qp_in, qp_out);

    // if (acados_return != ACADOS_SUCCESS)
    //     return -1;

    printf("STATUS: %d\n", acados_return);

//    print_dense_qp_out(qp_out);

    /************************************************
     * compute inf norm of residuals
     ************************************************/

    double res[4];
    dense_qp_res *qp_res = dense_qp_res_create(&dims);
    dense_qp_res_ws *res_ws = dense_qp_res_workspace_create(&dims);
    d_compute_res_dense_qp(qp_in, qp_out, qp_res, res_ws);

#if 0
	blasfeo_print_exp_tran_dvec(dims.nv+2*dims.ns, qp_res->res_g, 0);
	blasfeo_print_exp_tran_dvec(dims.ne, qp_res->res_b, 0);
	blasfeo_print_exp_tran_dvec(2*dims.nb+2*dims.ng+2*dims.ns, qp_res->res_d, 0);
	blasfeo_print_exp_tran_dvec(2*dims.nb+2*dims.ng+2*dims.ns, qp_res->res_m, 0);
#endif

    dense_qp_inf_norm_residuals(&dims, qp_in, qp_out, res);
    printf("\ninf norm res: %e, %e, %e, %e\n\n", res[0], res[1], res[2], res[3]);

    free(config);
    free(opts);
    free(qp_in);
    free(qp_out);
    free(qp_solver);
}

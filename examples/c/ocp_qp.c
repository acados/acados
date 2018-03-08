#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/print.h"
#include "acados_c/ocp_qp_interface.h"

#define QP_HORIZON 5

int main() {

    double A[] = {1, 0, 1, 1};
    double B[] = {0, 1};
    double b[] = {0, 0};

    double Q[] = {1, 0, 0, 1};
    double S[] = {0, 0};
    double R[] = {1};
    double q[] = {1, 1};
    double r[] = {0};

    double x0[] = {1, 1};
    int idxb0[] = {1, 2};

    ocp_qp_solver_plan plan;
    plan.qp_solver = PARTIAL_CONDENSING_HPIPM;

    ocp_qp_xcond_solver_config *config = ocp_qp_config_create(&plan, QP_HORIZON);

    ocp_qp_dims *dims = ocp_qp_dims_create(QP_HORIZON);

    dims->N = QP_HORIZON;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        dims->nx[ii] = 2;
        dims->ng[ii] = 0;
        dims->ns[ii] = 0;
        dims->nbu[ii] = 0;

        if (ii == 0)
            dims->nbx[ii] = 2;
        else
            dims->nbx[ii] = 0;

        dims->nb[ii] =  dims->nbx[ii] +  dims->nbu[ii];

        if (ii < dims->N)
            dims->nu[ii] = 1;
        else
            dims->nu[ii] = 0;
    }

    ocp_qp_in *qp_in = ocp_qp_in_create(config, dims);

    double *hA[] = {A, A, A, A, A};
    double *hB[] = {B, B, B, B, B};
    double *hb[] = {b, b, b, b, b};
    double *hQ[] = {Q, Q, Q, Q, Q, Q};
    double *hS[] = {S, S, S, S, S, S};
    double *hR[] = {R, R, R, R, R};
    double *hq[] = {q, q, q, q, q, q};
    double *hr[] = {r, r, r, r, r, r};
    int *hidxb[] = {idxb0};
    double *hlb[] = {x0};
    double *hub[] = {x0};
    double *hC[] = {};
    double *hD[] = {};
    double *hlg[] = {};
    double *hug[] = {};

    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, qp_in);

    print_ocp_qp_in(qp_in);

    void *opts = ocp_qp_opts_create(config, dims);

    ocp_qp_out *qp_out = ocp_qp_out_create(config, dims);

    // TODO(dimitris): only have N2 in one place!!
    // printf("N2 in config = %d\n", config->N2);
    // printf("N2 in opts = %d\n", ((ocp_qp_partial_condensing_args *)(((ocp_qp_partial_condensing_solver_opts *)opts)->pcond_opts))->N2);

    ocp_qp_solver *qp_solver = ocp_qp_create(config, dims, opts);

    int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

    if (acados_return != ACADOS_SUCCESS)
        return -1;

    print_ocp_qp_out(qp_out);

    ocp_qp_free(qp_solver, qp_in, qp_out);
}

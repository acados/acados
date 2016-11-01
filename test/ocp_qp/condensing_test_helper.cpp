#include "test/ocp_qp/condensing_test_helper.h"
#include "test/test_utils/eigen.h"
#include "test/test_utils/read_matrix.h"
#include "test/test_utils/zeros.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

void calculate_num_state_bounds(const ocp_qp_in *in, condensing_workspace *work) {
    i_zeros(&work->nstate_bounds, in->N+1, 1);
    int_t num_state_bounds;
    for (int_t i = 1; i <= in->N; i++) {
        num_state_bounds = 0;
        for (int_t j = 0; j < in->nb[i]; j++) {
            if (in->idxb[i][j] < in->nx[i])
                num_state_bounds++;
        }
        work->nstate_bounds[i] = num_state_bounds;
    }
}

int_t get_num_condensed_vars(const ocp_qp_in *in) {
    int_t num_condensed_vars = 0;
    // TODO(robin): this only holds for MPC, not MHE
    num_condensed_vars += 0*(in->nx[1]);
    for (int_t i = 0; i < in->N; i++)
        num_condensed_vars += in->nu[i];
    return num_condensed_vars;
}

int_t get_num_constraints(const ocp_qp_in *in, condensing_workspace *work) {
    calculate_num_state_bounds(in, work);
    int_t num_constraints = in->nc[0];
    for (int_t i = 1; i <= in->N; i++) {
        num_constraints += in->nc[i] + work->nstate_bounds[i];
    }
    return num_constraints;
}

void allocateForUnconstrainedQPData(int_t N, int_t nx, int_t nu, ocp_qp_in * const qp) {
    qp->N = N;
    i_zeros((int_t **) &qp->nx, N+1, 1);
    i_zeros((int_t **) &qp->nu, N, 1);
    i_zeros((int_t **) &qp->nb, N+1, 1);
    i_zeros((int_t **) &qp->nc, N+1, 1);
    qp->lb = (const real_t **) malloc(sizeof(*qp->lb));
    d_zeros((real_t **) &qp->lb[0], nx, 1);
    qp->ub = (const real_t **) malloc(sizeof(*qp->ub));
    d_zeros((real_t **) &qp->ub[0], nx, 1);
    qp->A = (const real_t **) malloc(sizeof(*qp->A) * N);
    qp->B = (const real_t **) malloc(sizeof(*qp->B) * N);
    qp->b = (const real_t **) malloc(sizeof(*qp->b) * N);
    qp->Q = (const real_t **) malloc(sizeof(*qp->Q) * (N+1));
    qp->S = (const real_t **) malloc(sizeof(*qp->S) * N);
    qp->R = (const real_t **) malloc(sizeof(*qp->R) * N);
    qp->q = (const real_t **) malloc(sizeof(*qp->q) * (N+1));
    qp->r = (const real_t **) malloc(sizeof(*qp->r) * N);
    for (int_t i = 0; i < N; i++) {
        d_zeros((real_t **) &qp->A[i], nx, nx);
        d_zeros((real_t **) &qp->B[i], nx, nu);
        d_zeros((real_t **) &qp->b[i], nx, 1);
        d_zeros((real_t **) &qp->Q[i], nx, nx);
        d_zeros((real_t **) &qp->S[i], nu, nx);
        d_zeros((real_t **) &qp->R[i], nu, nu);
        d_zeros((real_t **) &qp->q[i], nx, 1);
        d_zeros((real_t **) &qp->r[i], nu, 1);
    }
    d_zeros((real_t **) &qp->Q[N], nx, nx);
    d_zeros((real_t **) &qp->q[N], nx, 1);
}

void allocateForCondensingData(const ocp_qp_in * const qp_in, condensing_in *in,
    condensing_out *out, condensing_workspace *work) {

    int_t N = qp_in->N;
    const int_t *nc = qp_in->nc;

    // condensing input
    in->qp_input = (ocp_qp_in *) qp_in;

    // condensing output
    int_t nconvars = get_num_condensed_vars(qp_in);
    int_t nconstraints = get_num_constraints(qp_in, work);
    d_zeros(&out->H, nconvars, nconvars);
    d_zeros(&out->h, nconvars, 1);
    d_zeros(&out->lb, nconvars, 1);
    d_zeros(&out->ub, nconvars, 1);
    d_zeros(&out->A, nconstraints, nconvars);
    d_zeros(&out->lbA, nconstraints, 1);
    d_zeros(&out->ubA, nconstraints, 1);

    real_t QPOASES_INFTY = 1.0e8;
    for (int_t i = 0; i < nconvars; i++) {
        out->lb[i] = -QPOASES_INFTY;
        out->ub[i] = +QPOASES_INFTY;
    }
    for (int_t i = 0; i < nconstraints; i++) {
        out->lbA[i] = -QPOASES_INFTY;
        out->ubA[i] = +QPOASES_INFTY;
    }

    // condensing workspace
    work->nconvars = nconvars;
    work->nconstraints = nconstraints;
    work->G = (real_t ***) malloc(sizeof(*work->G) * N);
    work->g = (real_t **) malloc(sizeof(*work->g) * N);
    work->D = (real_t ***) malloc(sizeof(*work->D) * (N+1));
    for (int_t i = 0; i < N; i++) {
        work->G[i] = (real_t **) malloc(sizeof(*(work->G[i])) * (i+1));
        work->D[i] = (real_t **) malloc(sizeof(*(work->D[i])) * (i+1));
        d_zeros(&work->g[i], qp_in->nx[i], 1);
        for (int_t j = 0; j <= i; j++) {
            d_zeros(&work->G[i][j], qp_in->nx[i], qp_in->nu[j]);
            d_zeros(&work->D[i][j], nc[i], qp_in->nu[j]);
        }
    }
    work->D[N] = (real_t **) malloc(sizeof(*(work->D[N])) * N);
    for (int_t i = 0; i < N; i++) {
        d_zeros(&work->D[N][i], nc[N], qp_in->nu[i]);
    }
    d_zeros(&work->W1_x, qp_in->nx[0], qp_in->nx[0]);
    d_zeros(&work->W2_x, qp_in->nx[0], qp_in->nx[0]);
    d_zeros(&work->W1_u, qp_in->nx[0], qp_in->nu[0]);
    d_zeros(&work->W2_u, qp_in->nx[0], qp_in->nu[0]);
    d_zeros(&work->w1, qp_in->nx[0], 1);
    d_zeros(&work->w2, qp_in->nx[0], 1);
    d_zeros(&work->Sx0, qp_in->nu[0], 1);
}

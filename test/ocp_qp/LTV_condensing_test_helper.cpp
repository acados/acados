#include "test/ocp_qp/LTI_condensing_test_helper.h"
#include "test/ocp_qp/condensing_test_helper.h"
#include "test/test_utils/read_matrix.h"
#include "test/test_utils/zeros.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

extern void readLTVInputDimensionsFromFile(int_t *N, int_t *nx, int_t *nu);

void readTVUnconstrainedInputDataFromFile(int_t N, int_t nx, int_t nu, MatrixXd *A, MatrixXd *B,
    MatrixXd *b, VectorXd *x0, MatrixXd *Q, MatrixXd *S, MatrixXd *R, MatrixXd *q, MatrixXd *r) {

    *A = readMatrixFromFile("LTV/A.dat", nx, N*nx);
    *B = readMatrixFromFile("LTV/B.dat", nx, N*nu);
    *b = readMatrixFromFile("LTV/bv.dat", nx, N);
    *x0 = readVectorFromFile("LTV/x0.dat", nx);
    *Q = readMatrixFromFile("LTV/Q.dat", nx, (N+1)*nx);
    *S = readMatrixFromFile("LTV/S.dat", nu, N*nx);
    *R = readMatrixFromFile("LTV/R.dat", nu, N*nu);
    *q = readMatrixFromFile("LTV/qv.dat", nx, N+1);
    *r = readMatrixFromFile("LTV/rv.dat", nu, N);
}

void fillWithTVUnconstrainedData(ocp_qp_in *qp, VectorXd *x0) {
    int_t N, nx, nu;
    MatrixXd A, B, b, Q, S, R, q, r;

    readLTVInputDimensionsFromFile(&N, &nx, &nu);
    readTVUnconstrainedInputDataFromFile(N, nx, nu, &A, &B, &b, x0, &Q, &S, &R, &q, &r);
    allocateForUnconstrainedQPData(N, nx, nu, qp);

    int_t nx_vector[N+1], nu_vector[N], nb_vector[N+1], nc_vector[N+1];
    // Initial state
    nb_vector[0] = nx;
    for (int_t i = 0; i < N; i++) {
        nx_vector[i] = nx;
        nu_vector[i] = nu;
        nb_vector[i] = 0;
        nc_vector[i] = 0;
        memcpy((void *) qp->A[i], (void *) A.block(0, nx*i, nx, nx).data(),
            sizeof(*qp->A[i]) * (nx*nx));
        memcpy((void *) qp->B[i], (void *) B.block(0, nu*i, nx, nu).data(),
            sizeof(*qp->B[i]) * (nx*nu));
        memcpy((void *) qp->b[i], (void *) b.block(0, i, nx, 1).data(),
            sizeof(*qp->b[i]) * nx);
        memcpy((void *) qp->Q[i], (void *) Q.block(0, nx*i, nx, nx).data(),
            sizeof(*qp->Q[i]) * (nx*nx));
        memcpy((void *) qp->S[i], (void *) S.block(0, nx*i, nu, nx).data(),
            sizeof(*qp->S[i]) * (nu*nx));
        memcpy((void *) qp->R[i], (void *) R.block(0, nu*i, nu, nu).data(),
            sizeof(*qp->R[i]) * (nu*nu));
        memcpy((void *) qp->q[i], (void *) q.block(0, i, nx, 1).data(),
            sizeof(*qp->q[i]) * nx);
        memcpy((void *) qp->r[i], (void *) r.block(0, i, nu, 1).data(),
            sizeof(*qp->r[i]) * nu);
    }
    // Final state
    nx_vector[N] = nx;
    nb_vector[N] = 0;
    nc_vector[N] = 0;
    memcpy((void *) qp->Q[N], (void *) Q.block(0, nx*N, nx, nx).data(),
        sizeof(*qp->Q[N]) * (nx*nx));
    memcpy((void *) qp->q[N], (void *) q.block(0, N, nx, 1).data(),
        sizeof(*qp->q[N]) * nx);

    // Copy the dimensions
    memcpy((void *) qp->nx, (void *) &nx_vector[0], sizeof(*qp->nx) * (N+1));
    memcpy((void *) qp->nu, (void *) &nu_vector[0], sizeof(*qp->nu) * N);
    memcpy((void *) qp->nb, (void *) &nb_vector[0], sizeof(*qp->nb) * (N+1));
    memcpy((void *) qp->nc, (void *) &nc_vector[0], sizeof(*qp->nc) * (N+1));
    // Copy the initial constraint
    memcpy((void *) qp->lb[0], (void *) (*x0).data(), sizeof(*qp->lb[0]) * nx);
    memcpy((void *) qp->ub[0], (void *) (*x0).data(), sizeof(*qp->ub[0]) * nx);
}

void readTVBoundsDataFromFile(int_t N, int_t nx, int_t nu, MatrixXd *lbwx,
    MatrixXd *ubwx, MatrixXd *lbwu, MatrixXd *ubwu) {
    *lbwx = readMatrixFromFile("LTV/lower_bound_x.dat", nx, N+1);
    *ubwx = readMatrixFromFile("LTV/upper_bound_x.dat", nx, N+1);
    *lbwu = readMatrixFromFile("LTV/lower_bound_u.dat", nu, N);
    *ubwu = readMatrixFromFile("LTV/upper_bound_u.dat", nu, N);
}

void fillWithTVBoundsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu) {
    MatrixXd lbwx, ubwx, lbwu, ubwu;
    int_t nb_vector[N+1];
    int_t nb = (int_t) readMatrix("LTV/nb.dat")(0, 0);
    VectorXd lbw_block(nb);
    VectorXd ubw_block(nb);
    Eigen::VectorXi idxb_vector(nb);

    readTVBoundsDataFromFile(N, nx, nu, &lbwx, &ubwx, &lbwu, &ubwu);
    qp->lb = (const real_t **) malloc(sizeof(*qp->lb) * (N+1));
    qp->ub = (const real_t **) malloc(sizeof(*qp->ub) * (N+1));
    qp->idxb = (const int_t **) malloc(sizeof(*qp->idxb) * (N+1));
    for (int_t i = 0; i < nb; i++) {
        idxb_vector(i) = i;
    }
    for (int_t i = 0; i < N; i++) {
        nb_vector[i] = nb;
        lbw_block << lbwx.block(0, i, nx, 1), lbwu.block(0, i, nu, 1);
        ubw_block << ubwx.block(0, i, nx, 1), ubwu.block(0, i, nu, 1);
        i_zeros((int_t **) &qp->idxb[i], nb_vector[i], 1);
        memcpy((void *) qp->idxb[i], (void *) idxb_vector.data(),
            sizeof(*qp->idxb[i]) * nb_vector[i]);
        d_zeros((real_t **) &qp->lb[i], nb, 1);
        memcpy((void *) qp->lb[i], (void *) lbw_block.data(),
            sizeof(*qp->lb[i]) * (nb));
        d_zeros((real_t **) &qp->ub[i], nb, 1);
        memcpy((void *) qp->ub[i], (void *) ubw_block.data(),
            sizeof(*qp->ub[i]) * (nb));
    }
    nb_vector[N] = nx;
    i_zeros((int_t **) &qp->idxb[N], nb_vector[N], 1);
    memcpy((void *) qp->idxb[N], (void *) idxb_vector.data(),
        sizeof(*qp->idxb[N]) * nb_vector[N]);
    memcpy((void *) qp->nb, (void *) nb_vector,
        sizeof(*qp->nb) * (N+1));
    d_zeros((real_t **) &qp->lb[N], nb, 1);
    memcpy((void *) qp->lb[N], (void *) lbwx.block(0, N, nx, 1).data(),
        sizeof(*qp->lb[N]) * nx);
    d_zeros((real_t **) &qp->ub[N], nb, 1);
    memcpy((void *) qp->ub[N], (void *) ubwx.block(0, N, nx, 1).data(),
        sizeof(*qp->ub[N]) * nx);
}

void readTVGeneralConstraintsDataFromFile(int_t N, int_t nx, int_t nu, int_t nc,
    MatrixXd *Cx, MatrixXd *Cu, MatrixXd *lbc, MatrixXd *ubc) {

    *Cx = readMatrixFromFile("LTV/general_constraint_x.dat", nc, (N+1)*nx);
    *Cu = readMatrixFromFile("LTV/general_constraint_u.dat", nc, N*nu);
    *lbc = readMatrixFromFile("LTV/general_constraint_lb.dat", nc, N+1);
    *ubc = readMatrixFromFile("LTV/general_constraint_ub.dat", nc, N+1);
}

void fillWithTVGeneralConstraintsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu) {
    MatrixXd Cx, Cu;
    MatrixXd lbc, ubc;
    int_t nc_vector[N+1];
    int_t nc = (int_t) readMatrix("LTV/nc.dat")(0, 0);

    readTVGeneralConstraintsDataFromFile(N, nx, nu, nc, &Cx, &Cu, &lbc, &ubc);
    qp->Cx = (const real_t **) malloc(sizeof(*qp->Cx) * (N+1));
    qp->Cu = (const real_t **) malloc(sizeof(*qp->Cu) * N);
    qp->lc = (const real_t **) malloc(sizeof(*qp->lc) * (N+1));
    qp->uc = (const real_t **) malloc(sizeof(*qp->uc) * (N+1));
    for (int_t i = 0; i < N; i++) {
        nc_vector[i] = nc;
        d_zeros((real_t **) &qp->Cx[i], nc, nx);
        memcpy((void *) qp->Cx[i], (void *) Cx.block(0, i*nx, nc, nx).data(),
            sizeof(*qp->Cx[i]) * (nc*nx));
        d_zeros((real_t **) &qp->Cu[i], nc, nu);
        memcpy((void *) qp->Cu[i], (void *) Cu.block(0, i*nu, nc, nu).data(),
            sizeof(*qp->Cu[i]) * (nc*nu));
        d_zeros((real_t **) &qp->lc[i], nc, 1);
        memcpy((void *) qp->lc[i], (void *) lbc.block(0, i, nc, 1).data(),
            sizeof(*qp->lc[i]) * nc);
        d_zeros((real_t **) &qp->uc[i], nc, 1);
        memcpy((void *) qp->uc[i], (void *) ubc.block(0, i, nc, 1).data(),
            sizeof(*qp->uc[i]) * nc);
    }
    nc_vector[N] = nc;
    d_zeros((real_t **) &qp->Cx[N], nc, nx);
    memcpy((void *) qp->Cx[N], (void *) Cx.block(0, N*nx, nc, nx).data(),
        sizeof(*qp->Cx[N]) * (nc*nx));
    d_zeros((real_t **) &qp->lc[N], nc, 1);
    memcpy((void *) qp->lc[N], (void *) lbc.block(0, N, nc, 1).data(),
        sizeof(*qp->lc[N]) * nc);
    d_zeros((real_t **) &qp->uc[N], nc, 1);
    memcpy((void *) qp->uc[N], (void *) ubc.block(0, N, nc, 1).data(),
        sizeof(*qp->uc[N]) * nc);
    memcpy((void *) qp->nc, (void *) nc_vector, sizeof(*qp->nc) * N+1);
}

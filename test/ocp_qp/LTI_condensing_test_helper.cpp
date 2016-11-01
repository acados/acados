#include "test/ocp_qp/LTI_condensing_test_helper.h"
#include "test/ocp_qp/condensing_test_helper.h"
#include "test/test_utils/read_matrix.h"
#include "test/test_utils/zeros.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

extern void readLTIInputDimensionsFromFile(int_t *N, int_t *nx, int_t *nu);

void readUnconstrainedInputDataFromFile(int_t nx, int_t nu, MatrixXd *A, MatrixXd *B,
    VectorXd *b, VectorXd *x0, MatrixXd *Q, MatrixXd *S, MatrixXd *R, VectorXd *q, VectorXd *r) {

    *A = readMatrixFromFile("LTI/A.dat", nx, nx);
    *B = readMatrixFromFile("LTI/B.dat", nx, nu);
    *b = readVectorFromFile("LTI/bv.dat", nx);
    *x0 = readVectorFromFile("LTI/x0.dat", nx);
    *Q = readMatrixFromFile("LTI/Q.dat", nx, nx);
    *S = readMatrixFromFile("LTI/S.dat", nu, nx);
    *R = readMatrixFromFile("LTI/R.dat", nu, nu);
    *q = readVectorFromFile("LTI/qv.dat", nx);
    *r = readVectorFromFile("LTI/rv.dat", nu);
}

void fillWithUnconstrainedData(ocp_qp_in *qp, VectorXd *x0) {
    int_t N, nx, nu;
    MatrixXd A, B, Q, S, R;
    VectorXd b, q, r;

    readLTIInputDimensionsFromFile(&N, &nx, &nu);
    readUnconstrainedInputDataFromFile(nx, nu, &A, &B, &b, x0, &Q, &S, &R, &q, &r);
    allocateForUnconstrainedQPData(N, nx, nu, qp);

    int_t nx_vector[N+1], nu_vector[N], nb_vector[N+1], nc_vector[N+1];
    // Initial state
    nb_vector[0] = nx;
    for (int_t i = 0; i < N; i++) {
        nx_vector[i] = nx;
        nu_vector[i] = nu;
        nb_vector[i] = 0;
        nc_vector[i] = 0;
        memcpy((void *) qp->A[i], (void *) A.data(), sizeof(*qp->A[i]) * (nx*nx));
        memcpy((void *) qp->B[i], (void *) B.data(), sizeof(*qp->B[i]) * (nx*nu));
        memcpy((void *) qp->b[i], (void *) b.data(), sizeof(*qp->b[i]) * nx);
        memcpy((void *) qp->Q[i], (void *) Q.data(), sizeof(*qp->Q[i]) * (nx*nx));
        memcpy((void *) qp->S[i], (void *) S.data(), sizeof(*qp->S[i]) * (nu*nx));
        memcpy((void *) qp->R[i], (void *) R.data(), sizeof(*qp->R[i]) * (nu*nu));
        memcpy((void *) qp->q[i], (void *) q.data(), sizeof(*qp->q[i]) * nx);
        memcpy((void *) qp->r[i], (void *) r.data(), sizeof(*qp->r[i]) * nu);
    }
    // Final state
    nx_vector[N] = nx;
    nb_vector[N] = 0;
    nc_vector[N] = 0;
    memcpy((void *) qp->Q[N], (void *) Q.data(), sizeof(*qp->Q[N]) * (nx*nx));
    memcpy((void *) qp->q[N], (void *) q.data(), sizeof(*qp->q[N]) * nx);

    // Copy the dimensions
    memcpy((void *) qp->nx, (void *) &nx_vector[0], sizeof(*qp->nx) * (N+1));
    memcpy((void *) qp->nu, (void *) &nu_vector[0], sizeof(*qp->nu) * N);
    memcpy((void *) qp->nb, (void *) &nb_vector[0], sizeof(*qp->nb) * (N+1));
    memcpy((void *) qp->nc, (void *) &nc_vector[0], sizeof(*qp->nc) * (N+1));
    // Copy the initial constraint
    memcpy((void *) qp->lb[0], (void *) (*x0).data(), sizeof(*qp->lb[0]) * nx);
    memcpy((void *) qp->ub[0], (void *) (*x0).data(), sizeof(*qp->ub[0]) * nx);
}

void readBoundsDataFromFile(int_t nx, int_t nu, VectorXd *lbwx,
    VectorXd *ubwx, VectorXd *lbwu, VectorXd *ubwu) {
    *lbwx = readVectorFromFile("LTI/lower_bound_x.dat", nx);
    *ubwx = readVectorFromFile("LTI/upper_bound_x.dat", nx);
    *lbwu = readVectorFromFile("LTI/lower_bound_u.dat", nu);
    *ubwu = readVectorFromFile("LTI/upper_bound_u.dat", nu);
}

void fillWithBoundsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu) {
    VectorXd lbwx, ubwx, lbwu, ubwu;
    int_t nb_vector[N+1];
    int_t nb = (int_t) readMatrix("LTI/nb.dat")(0, 0);
    VectorXd lbw_block(nb);
    VectorXd ubw_block(nb);
    Eigen::VectorXi idxb_vector(nb);

    readBoundsDataFromFile(nx, nu, &lbwx, &ubwx, &lbwu, &ubwu);
    lbw_block << lbwx, lbwu;
    ubw_block << ubwx, ubwu;
    qp->lb = (const real_t **) malloc(sizeof(*qp->lb) * (N+1));
    qp->ub = (const real_t **) malloc(sizeof(*qp->ub) * (N+1));
    qp->idxb = (const int_t **) malloc(sizeof(*qp->idxb) * (N+1));
    for (int_t i = 0; i < nb; i++) {
        idxb_vector(i) = i;
    }
    for (int_t i = 0; i < N; i++) {
        nb_vector[i] = nb;
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
    memcpy((void *) qp->idxb[N], (void *) idxb_vector.data(), sizeof(*qp->idxb[N]) * nb_vector[N]);
    memcpy((void *) qp->nb, (void *) nb_vector, sizeof(*qp->nb) * (N+1));
    d_zeros((real_t **) &qp->lb[N], nb, 1);
    memcpy((void *) qp->lb[N], (void *) lbwx.data(), sizeof(*qp->lb[N]) * nx);
    d_zeros((real_t **) &qp->ub[N], nb, 1);
    memcpy((void *) qp->ub[N], (void *) ubwx.data(), sizeof(*qp->ub[N]) * nx);
}

void readGeneralConstraintsDataFromFile(int_t nx, int_t nu, int_t nc,
    MatrixXd *Cx, MatrixXd *Cu, VectorXd *lbc, VectorXd *ubc) {

    *Cx = readMatrixFromFile("LTI/general_constraint_x.dat", nc, nx);
    *Cu = readMatrixFromFile("LTI/general_constraint_u.dat", nc, nu);
    *lbc = readVectorFromFile("LTI/general_constraint_lb.dat", nc);
    *ubc = readVectorFromFile("LTI/general_constraint_ub.dat", nc);
}

void fillWithGeneralConstraintsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu) {
    MatrixXd Cx, Cu;
    VectorXd lbc, ubc;
    int_t nc_vector[N+1];
    int_t nc = (int_t) readMatrix("LTI/nc.dat")(0, 0);

    readGeneralConstraintsDataFromFile(nx, nu, nc, &Cx, &Cu, &lbc, &ubc);
    qp->Cx = (const real_t **) malloc(sizeof(*qp->Cx) * (N+1));
    qp->Cu = (const real_t **) malloc(sizeof(*qp->Cu) * N);
    qp->lc = (const real_t **) malloc(sizeof(*qp->lc) * (N+1));
    qp->uc = (const real_t **) malloc(sizeof(*qp->uc) * (N+1));
    for (int_t i = 0; i < N; i++) {
        nc_vector[i] = nc;
        d_zeros((real_t **) &qp->Cx[i], nc, nx);
        memcpy((void *) qp->Cx[i], (void *) Cx.data(), sizeof(*qp->Cx[i]) * (nc*nx));
        d_zeros((real_t **) &qp->Cu[i], nc, nu);
        memcpy((void *) qp->Cu[i], (void *) Cu.data(), sizeof(*qp->Cu[i]) * (nc*nu));
        d_zeros((real_t **) &qp->lc[i], nc, 1);
        memcpy((void *) qp->lc[i], (void *) lbc.data(), sizeof(*qp->lc[i]) * nc);
        d_zeros((real_t **) &qp->uc[i], nc, 1);
        memcpy((void *) qp->uc[i], (void *) ubc.data(), sizeof(*qp->uc[i]) * nc);
    }
    nc_vector[N] = nc;
    d_zeros((real_t **) &qp->Cx[N], nc, nx);
    memcpy((void *) qp->Cx[N], (void *) Cx.data(), sizeof(*qp->Cx[N]) * (nc*nx));
    d_zeros((real_t **) &qp->lc[N], nc, 1);
    memcpy((void *) qp->lc[N], (void *) lbc.data(), sizeof(*qp->lc[N]) * nc);
    d_zeros((real_t **) &qp->uc[N], nc, 1);
    memcpy((void *) qp->uc[N], (void *) ubc.data(), sizeof(*qp->uc[N]) * nc);
    memcpy((void *) qp->nc, (void *) nc_vector, sizeof(*qp->nc) * N+1);
}

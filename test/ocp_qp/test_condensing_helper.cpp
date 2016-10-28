#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

static void readInputDimensionsFromFile(int_t *N, int_t *nx, int_t *nu);
static MatrixXd readMatrixFromFile(std::string filename, int_t rows, int_t cols);
static VectorXd readVectorFromFile(std::string filename, int_t length);

static void i_zeros(int_t **pA, int_t row, int_t col) {
    void *temp = malloc((row*col)*sizeof(int_t));
    *pA = (int_t *) temp;
    int_t *A = *pA;
    int_t i;
    for (i = 0; i < row*col; i++) A[i] = 0;
}

static void d_zeros(real_t **pA, int_t row, int_t col) {
    void *temp = malloc((row*col)*sizeof(real_t));
    *pA = (real_t *) temp;
    real_t *A = *pA;
    int_t i;
    for (i = 0; i < row*col; i++) A[i] = 0.0;
}

static void calculate_num_state_bounds(ocp_qp_in const *in, condensing_workspace *work) {
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

static int_t get_num_condensed_vars(ocp_qp_in const *in) {
    int_t num_condensed_vars = 0;
    // TODO(robin): this only holds for MPC, not MHE
    num_condensed_vars += 0*(in->nx[1]);
    for (int_t i = 0; i < in->N; i++)
        num_condensed_vars += in->nu[i];
    return num_condensed_vars;
}

static int_t get_num_constraints(ocp_qp_in const *in, condensing_workspace *work) {
    calculate_num_state_bounds(in, work);
    int_t num_constraints = in->nc[0];
    for (int_t i = 1; i <= in->N; i++) {
        num_constraints += in->nc[i] + work->nstate_bounds[i];
    }
    return num_constraints;
}

static void fill_in_condensing_structs(ocp_qp_in const *qp_in, condensing_in *in,
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

static void readUnconstrainedInputDataFromFile(int_t nx, int_t nu, MatrixXd *A, MatrixXd *B,
    VectorXd *b, VectorXd *x0, MatrixXd *Q, MatrixXd *S, MatrixXd *R, VectorXd *q, VectorXd *r) {
    *A = readMatrixFromFile("A.dat", nx, nx);
    *B = readMatrixFromFile("B.dat", nx, nu);
    *b = readVectorFromFile("bv.dat", nx);
    *x0 = readVectorFromFile("x0.dat", nx);
    *Q = readMatrixFromFile("Q.dat", nx, nx);
    *S = readMatrixFromFile("S.dat", nu, nx);
    *R = readMatrixFromFile("R.dat", nu, nu);
    *q = readVectorFromFile("qv.dat", nx);
    *r = readVectorFromFile("rv.dat", nu);
}

// static void readBoundsInputDataFromFile(int_t N, int_t nx, int_t nu, VectorXd *lb, VectorXd *ub) {
//     *lb = readVectorFromFile("lower_bound.dat", N*(nx+nu)+nx);
//     *ub = readVectorFromFile("upper_bound.dat", N*(nx+nu)+nx);
// }

static void fillWithUnconstrainedData(ocp_qp_in *qp, condensing_in *input,
    condensing_out *output, condensing_workspace *work, VectorXd *x0) {
    int_t N, nx, nu;
    readInputDimensionsFromFile(&N, &nx, &nu);
    MatrixXd A, B, Q, S, R;
    VectorXd b, q, r;
    readUnconstrainedInputDataFromFile(nx, nu, &A, &B, &b, x0, &Q, &S, &R, &q, &r);

    int_t nx_vector[N+1], nu_vector[N], nb_vector[N+1], nc_vector[N+1];

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

    fill_in_condensing_structs(qp, input, output, work);
}

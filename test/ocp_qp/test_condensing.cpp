#include "catch/include/catch.hpp"

#include "test/test_utils/read_matrix.hpp"
#include "acados/ocp_qp/condensing.c"

#define EPS 1.0e-15

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

static void calculate_num_state_bounds(ocp_qp_in *in, condensing_workspace *work) {
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

static int_t get_num_condensed_vars(ocp_qp_in *in) {
    int_t num_condensed_vars = 0;
    // TODO(robin): this only holds for MPC, not MHE
    num_condensed_vars += 0*(in->nx[1]);
    for (int_t i = 0; i < in->N; i++)
        num_condensed_vars += in->nu[i];
    return num_condensed_vars;
}

static int_t get_num_constraints(ocp_qp_in *in, condensing_workspace *work) {
    calculate_num_state_bounds(in, work);
    int_t num_constraints = in->nc[0];
    for (int_t i = 1; i <= in->N; i++) {
        num_constraints += in->nc[i] + work->nstate_bounds[i];
    }
    return num_constraints;
}

static void fill_in_condensing_structs(ocp_qp_in *qp_in, condensing_in *in,
    condensing_out *out, condensing_workspace *work) {
    int_t N = qp_in->N;
    const int_t *nc = qp_in->nc;

    // condensing input
    in->qp_input = qp_in;

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

TEST_CASE("Unconstrained LTV system", "[condensing]") {
    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;

    // Read input and result from file
    int_t N = (int_t) readMatrix("N.dat")(0, 0);
    int_t nx = (int_t) readMatrix("nx.dat")(0, 0);
    int_t nu = (int_t) readMatrix("nu.dat")(0, 0);
    MatrixXd A = readMatrix("A.dat");
    REQUIRE(A.rows() == nx);
    REQUIRE(A.cols() == nx);
    MatrixXd B = readMatrix("B.dat");
    REQUIRE(B.rows() == nx);
    REQUIRE(B.cols() == nu);
    VectorXd b = readMatrix("c.dat");
    REQUIRE(b.rows() == nx);
    VectorXd x0 = readMatrix("x0.dat");
    REQUIRE(x0.rows() == nx);

    // Fill in data
    int_t nx_vector[N+1];
    int_t nu_vector[N];
    int_t nb_vector[N+1];
    int_t nc_vector[N+1];
    real_t *A_vector[N];
    real_t *B_vector[N];
    real_t *b_vector[N];
    real_t *lb_vector[1];
    real_t *ub_vector[1];

    for (int_t i = 0; i < N; i++) {
        nx_vector[i] = nx;
        nu_vector[i] = nu;
        nb_vector[i] = 0;
        nc_vector[i] = 0;
        A_vector[i] = A.data();
        B_vector[i] = B.data();
        b_vector[i] = b.data();
    }
    // Initial state
    nb_vector[0] = nx;
    lb_vector[0] = x0.data();
    ub_vector[0] = x0.data();
    // Final state
    nx_vector[N] = nx;
    nb_vector[N] = 0;
    nc_vector[N] = 0;

    // Finalize the data
    qp.N = N;
    qp.nx = (const int_t *) nx_vector;
    qp.nu = (const int_t *) nu_vector;
    qp.nb = (const int_t *) nb_vector;
    qp.nc = (const int_t *) nc_vector;
    qp.A = (const real_t **) A_vector;
    qp.B = (const real_t **) B_vector;
    qp.b = (const real_t **) b_vector;
    qp.lb = (const real_t **) lb_vector;
    qp.ub = (const real_t **) ub_vector;

    fill_in_condensing_structs(&qp, &input, &output, &work);

    SECTION("Transition vector") {
        calculate_transition_vector(&qp, &work, x0.data());
        VectorXd true_g = readMatrix("transition_vector.dat");
        VectorXd acados_g = VectorXd(true_g).setZero();
        for (int_t i = 0; i < N; i++) {
            acados_g.block(i*nx, 0, nx, 1) = Eigen::Map<MatrixXd>(work.g[i], nx, 1);
        }
        REQUIRE(acados_g.isApprox(true_g, EPS));
    }

    SECTION("Transition matrix") {
        calculate_transition_matrix(&qp, &work);
        MatrixXd true_G = readMatrix("transition_matrix.dat");
        MatrixXd acados_G = MatrixXd(true_G).setZero();
        for (int_t i = 0; i < N; i++) {
            for (int_t j = 0; j <= i; j++) {
                acados_G.block(i*nx, j*nu, nx, nu) = Eigen::Map<MatrixXd>(work.G[i][j], nx, nu);
            }
        }
        REQUIRE(acados_G.isApprox(true_G));
    }
}

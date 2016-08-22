// define number of repetitions
#define NREP 10

// system headers
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

// ACADOS headers
#include "acados/ocp_qp_condensing_qpoases.h"
#include "acados/tools.h"

#include "hpmpc/include/aux_d.h"

/************************************************
Mass-spring system: nx/2 masses connected each other with springs (in a row),
and the first and the last one to walls. nu (<=nx) controls act on the first nu
masses. The system is sampled with sampling time Ts.
************************************************/
void mass_spring_system(double Ts, int nx, int nu, double *A, double *B,
                        double *b, double *x0) {
    int nx2 = nx * nx;

    int info = 0;

    int pp = nx / 2;  // number of masses

    /************************************************
    * build the continuous time system
    ************************************************/

    double *T;
    d_zeros(&T, pp, pp);
    int ii;
    for (ii = 0; ii < pp; ii++) T[ii * (pp + 1)] = -2;
    for (ii = 0; ii < pp - 1; ii++) T[ii * (pp + 1) + 1] = 1;
    for (ii = 1; ii < pp; ii++) T[ii * (pp + 1) - 1] = 1;

    double *Z;
    d_zeros(&Z, pp, pp);
    double *I;
    d_zeros(&I, pp, pp);
    for (ii = 0; ii < pp; ii++) I[ii * (pp + 1)] = 1.0;  // = eye(pp);
    double *Ac;
    d_zeros(&Ac, nx, nx);
    dmcopy(pp, pp, Z, pp, Ac, nx);
    dmcopy(pp, pp, T, pp, Ac + pp, nx);
    dmcopy(pp, pp, I, pp, Ac + pp * nx, nx);
    dmcopy(pp, pp, Z, pp, Ac + pp * (nx + 1), nx);
    d_free(T);
    d_free(Z);
    d_free(I);

    d_zeros(&I, nu, nu);
    for (ii = 0; ii < nu; ii++) I[ii * (nu + 1)] = 1.0;  // I = eye(nu);
    double *Bc;
    d_zeros(&Bc, nx, nu);
    dmcopy(nu, nu, I, nu, Bc + pp, nx);
    d_free(I);

    /************************************************
    * compute the discrete time system
    ************************************************/

    double *bb;
    d_zeros(&bb, nx, 1);
    dmcopy(nx, 1, bb, nx, b, nx);

    dmcopy(nx, nx, Ac, nx, A, nx);
    dscal_3l(nx2, Ts, A);
    expm(nx, A);

    double *R, *J;
    d_zeros(&R, nx, nx);
    d_zeros(&J, nx, nx);
    for (ii = 0; ii < nx; ii++) J[ii * (nx + 1)] = 1.0;  // I = eye(nx);
    dmcopy(nx, nx, A, nx, R, nx);
    daxpy_3l(nx2, -1.0, J, R);
    dgemm_nn_3l(nx, nu, nx, R, nx, Bc, nx, B, nx);

    int *ipiv = (int *)malloc(nx * sizeof(int));
    dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
    free(ipiv);

    d_free(Ac);
    d_free(Bc);
    d_free(bb);
    d_free(R);
    d_free(J);

    /************************************************
    * initial state
    ************************************************/

    for (ii = 0; ii < nx; ii++) x0[ii] = 0;
    x0[0] = 2.5;
    x0[1] = 2.5;
}

int main() {

    int ii, jj;
    int nrep = NREP;

    int nx = 8;  // number of states (it has to be even for the mass-spring
                  // system test problem)
    int nu = 3;  // number of inputs (controllers) (it has to be at least 1 and
                  // at most nx/2 for the mass-spring system test problem)
    int N = 20;   // horizon length
    int nb = 11;  // number of box constrained inputs and states
    int ng = 0;  // 4;  // number of general constraints
    int ngN = 8;  // 4;  // number of general constraints at the last stage

    // int nbu = nu < nb ? nu : nb;
    int nbx = nb - nu > 0 ? nb - nu : 0;

    // stage-wise variant size
    int nxx[N + 1];
    for (ii = 0; ii <= N; ii++) nxx[ii] = nx;

    int nuu[N + 1];
    for (ii = 0; ii < N; ii++) nuu[ii] = nu;
    nuu[N] = 0;

    int nbb[N + 1];
    for (ii = 0; ii < N; ii++) nbb[ii] = nb;
    nbb[N] = nbx;

    int ngg[N + 1];
    for (ii = 0; ii < N; ii++) ngg[ii] = ng;
    ngg[0] = 1;
    ngg[N] = ngN;

    printf(
        " Test problem: mass-spring system with %d masses and %d controls.\n",
        nx / 2, nu);
    printf("\n");
    printf(
        " MPC problem size: %d states, %d inputs, %d horizon length, %d "
        "two-sided box constraints, %d two-sided general constraints.\n",
        nx, nu, N, nb, ng);
    printf("\n");

    /************************************************
    * dynamical system
    ************************************************/

    // state space matrices & initial state
    double *A;
    d_zeros(&A, nx, nx);  // states update matrix
    double *B;
    d_zeros(&B, nx, nu);  // inputs matrix
    double *b;
    d_zeros(&b, nx, 1);  // states offset
    double *x0;
    d_zeros(&x0, nx, 1);  // initial state

    // mass-spring system
    double Ts = 0.5;  // sampling time
    mass_spring_system(Ts, nx, nu, A, B, b, x0);

    for (jj = 0; jj < nx; jj++) b[jj] = 0.1;

    // compute b0 = b + A*x0
    double *b0;
    d_zeros(&b0, nx, 1);
    dcopy_3l(nx, b, 1, b0, 1);

    // then A0 is a matrix of size 0x0
    double *A0;
    d_zeros(&A0, 0, 0);

    /************************************************
    * box constraints
    ************************************************/

    int *idxb0;
    i_zeros(&idxb0, nbb[0], 1);
    double *lb0;
    d_zeros(&lb0, nbb[0], 1);
    double *ub0;
    d_zeros(&ub0, nbb[0], 1);
    for (jj = 0; jj < nx; jj++) {
        lb0[jj] = x0[jj];   //   xmin
        ub0[jj] = x0[jj];   //   xmax
        idxb0[jj] = jj;
    }
    for (; jj < nb; jj++) {
        lb0[jj] = -0.5;  //   umin
        ub0[jj] = +0.5;   //   umax
        idxb0[jj] = jj;
    }

    int *idxb1;
    i_zeros(&idxb1, nbb[1], 1);
    double *lb1;
    d_zeros(&lb1, nbb[1], 1);
    double *ub1;
    d_zeros(&ub1, nbb[1], 1);
    for (jj = 0; jj < nbx; jj++) {
        lb1[jj] = -4.0;  //   xmin
        ub1[jj] = +4.0;   //   xmax
        idxb1[jj] = jj;
    }
    for (; jj < nb; jj++) {
        lb1[jj] = -0.5;  //   umin
        ub1[jj] = +0.5;   //   umax
        idxb1[jj] = jj;
    }

    int *idxbN;
    i_zeros(&idxbN, nbb[N], 1);
    double *lbN;
    d_zeros(&lbN, nbb[N], 1);
    double *ubN;
    d_zeros(&ubN, nbb[N], 1);
    for (jj = 0; jj < nbx; jj++) {
        lbN[jj] = -4.0;  //   umin
        ubN[jj] = +4.0;   //   umax
        idxbN[jj] = jj;
    }

    /************************************************
    * general constraints
    ************************************************/

    double *C0;
    d_zeros(&C0, 1, nx);
    C0[0] = 1;
    double *D0;
    d_zeros(&D0, 1, nu);
    D0[0] = 1;
    double *lg0;
    d_zeros(&lg0, 1, 1);
    lg0[0] = 2.5;
    double *ug0;
    d_zeros(&ug0, 1, 1);
    ug0[0] = 5.5;

    double *C;
    d_zeros(&C, 0, nx);
    double *D;
    d_zeros(&D, 0, nu);
    double *lg;
    d_zeros(&lg, 0, 1);
    double *ug;
    d_zeros(&ug, 0, 1);

    double *CN;
    d_zeros(&CN, ngN, nx);
    for (ii = 0; ii < nx; ii++) CN[ii * (ngN + 1)] = 1.0;
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

    /************************************************
    * cost function
    ************************************************/

    double *Q;
    d_zeros(&Q, nx, nx);
    for (ii = 0; ii < nx; ii++) Q[ii * (nx + 1)] = 1.0;

    double *R;
    d_zeros(&R, nu, nu);
    for (ii = 0; ii < nu; ii++) R[ii * (nu + 1)] = 2.0;

    double *S;
    d_zeros(&S, nu, nx);

    double *q;
    d_zeros(&q, nx, 1);
    for (ii = 0; ii < nx; ii++) q[ii] = 0.1;

    double *r;
    d_zeros(&r, nu, 1);
    for (ii = 0; ii < nu; ii++) r[ii] = 0.2;

    /************************************************
    * problems data
    ************************************************/

    ocp_qp_input qp_in;
    ocp_qp_output qp_out;

    double *hA[N];
    double *hB[N];
    double *hb[N];
    double *hQ[N + 1];
    double *hS[N];
    double *hR[N];
    double *hq[N + 1];
    double *hr[N];
    double *hlb[N + 1];
    double *hub[N + 1];
    int *hidxb[N + 1];
    double *hC[N + 1];
    double *hD[N];
    double *hlg[N + 1];
    double *hug[N + 1];

    hA[0] = A;
    hB[0] = B;
    hb[0] = b;
    hQ[0] = Q;
    hS[0] = S;
    hR[0] = R;
    hq[0] = q;
    hr[0] = r;
    hlb[0] = lb0;
    hub[0] = ub0;
    hidxb[0] = idxb0;
    hC[0] = C0;
    hD[0] = D0;
    hlg[0] = lg0;
    hug[0] = ug0;
    for (ii = 1; ii < N; ii++) {
        hA[ii] = A;
        hB[ii] = B;
        hb[ii] = b;
        hQ[ii] = Q;
        hS[ii] = S;
        hR[ii] = R;
        hq[ii] = q;
        hr[ii] = r;
        hlb[ii] = lb1;
        hub[ii] = ub1;
        hidxb[ii] = idxb1;
        hC[ii] = C;
        hD[ii] = D;
        hlg[ii] = lg;
        hug[ii] = ug;
    }
    hQ[N] = Q;  // or maybe initialize to the solution of the DARE???
    hq[N] = q;  // or maybe initialize to the solution of the DARE???
    hlb[N] = lbN;
    hub[N] = ubN;
    hidxb[N] = idxbN;
    hC[N] = CN;
    hlg[N] = lgN;
    hug[N] = ugN;

    /************************************************
    * solution
    ************************************************/

    double *hx[N + 1];
    double *hu[N + 1];

    for (ii = 0; ii < N; ii++) {
        d_zeros(&hx[ii], nxx[ii], 1);
        d_zeros(&hu[ii], nuu[ii], 1);
    }
    d_zeros(&hx[N], nxx[N], 1);

    qp_in.N = N;
    qp_in.nx = nxx;
    qp_in.nu = nuu;
    qp_in.nb = nbb;
    qp_in.nc = ngg;
    qp_in.A = hA;
    qp_in.B = hB;
    qp_in.b = hb;
    qp_in.Q = hQ;
    qp_in.S = hS;
    qp_in.R = hR;
    qp_in.q = hq;
    qp_in.r = hr;
    qp_in.idxb = hidxb;
    qp_in.lb = hlb;
    qp_in.ub = hub;
    qp_in.Cx = hC;
    qp_in.Cu = hD;
    qp_in.lc = hlg;
    qp_in.uc = hug;

    qp_out.x = hx;
    qp_out.u = hu;

    /************************************************
    * solver arguments
    ************************************************/

    // solver arguments
    ocp_qp_condensing_qpoases_args args;
    args.dummy = 42.0;

    /************************************************
    * work space
    ************************************************/

    // int work_space_size = ocp_qp_condensing_qpoases_workspace_size(&qp_in, &args);
    // printf("\nwork space size: %d bytes\n", work_space_size);

    // double *work = (double *)malloc(work_space_size);

    /************************************************
    * call the solver
    ************************************************/

    int return_value;
    initialise_qpoases(&qp_in);

    struct timeval tv0, tv1;
    gettimeofday(&tv0, NULL);  // stop

    for (int rep = 0; rep < nrep; rep++) {
        // call the QP OCP solver
        return_value = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, NULL);
    }

    gettimeofday(&tv1, NULL);  // stop

    double time = (tv1.tv_sec - tv0.tv_sec) / (nrep + 0.0) +
                  (tv1.tv_usec - tv0.tv_usec) / (nrep * 1e6);

    printf("\nu = \n");
    for (ii = 0; ii < N; ii++) d_print_mat(1, nuu[ii], hu[ii], 1);

    printf("\nx = \n");
    for (ii = 0; ii <= N; ii++) d_print_mat(1, nxx[ii], hx[ii], 1);

    printf("\n");
    printf(" Average solution time over %d runs: %5.2e seconds\n", nrep, time);
    printf("\n\n");

    if (return_value != 0) {
        printf("\n qpOASES error ! No. %d\n", return_value);
    }

    /************************************************
    * free memory
    ************************************************/

    d_free(A);
    d_free(B);
    d_free(b);
    d_free(x0);
    d_free(A0);
    d_free(b0);
    d_free(C0);
    d_free(D0);
    d_free(lg0);
    d_free(ug0);
    d_free(Q);
    d_free(S);
    d_free(R);
    d_free(q);
    d_free(r);
    i_free(idxb0);
    d_free(lb0);
    d_free(ub0);
    i_free(idxb1);
    d_free(lb1);
    d_free(ub1);
    i_free(idxbN);
    d_free(lbN);
    d_free(ubN);
    d_free(C);
    d_free(D);
    d_free(lg);
    d_free(ug);
    d_free(CN);
    d_free(lgN);
    d_free(ugN);

    for (ii = 0; ii < N; ii++) {
        d_free(hx[ii]);
        d_free(hu[ii]);
    }
    d_free(hx[N]);

    return 0;
}

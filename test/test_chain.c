#include <stdio.h>
#include <stdlib.h>
#if defined(__APPLE__)
#include <mach/mach_time.h>
#else
#include <sys/stat.h>
#include <sys/time.h>
#endif
#include <math.h>
#include "hpmpc/include/aux_d.h"
#include "acados/utils/types.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "external/hpmpc/include/aux_d.h"
#include "test/casadi_chain/Chain_model.h"

#define NN 10
#define UMAX 2

#if defined(__APPLE__)
typedef struct acado_timer_ {
    uint64_t tic;
    uint64_t toc;
    mach_timebase_info_data_t tinfo;
} acado_timer;

static void acado_tic(acado_timer* t) {
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

static real_t acado_toc(acado_timer* t) {
    uint64_t duration; /* elapsed time in clock cycles*/

    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (real_t)duration / 1e9;
}
#else

typedef struct acado_timer_ {
    struct timeval tic;
    struct timeval toc;
} acado_timer;

static void acado_tic(acado_timer* t) {
    /* read current clock cycles */
    gettimeofday(&t->tic, 0);
}

static real_t acado_toc(acado_timer* t) {
    struct timeval temp;

    gettimeofday(&t->toc, 0);

    if ((t->toc.tv_usec - t->tic.tv_usec) < 0) {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
        temp.tv_usec = 1000000 + t->toc.tv_usec - t->tic.tv_usec;
    } else {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
        temp.tv_usec = t->toc.tv_usec - t->tic.tv_usec;
    }

    return (real_t)temp.tv_sec + (real_t)temp.tv_usec / 1e6;
}
#endif

// static void shift_states(real_t *w, real_t *x_end, int_t N) {
//    for (int_t i = 0; i < N; i++) {
//        for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = w[(i+1)*(NX+NU)+j];
//    }
//    for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] = x_end[j];
// }
//
// static void shift_controls(real_t *w, real_t *u_end, int_t N) {
//    for (int_t i = 0; i < N-1; i++) {
//        for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = w[(i+1)*(NX+NU)+NX+j];
//    }
//    for (int_t j = 0; j < NU; j++) w[(N-1)*(NX+NU)+NX+j] = u_end[j];
// }

// Simple SQP example for acados
int main() {
    int_t nil;

    for (int_t NMF = 1; NMF < 8; NMF++) {
        printf("\n------------ NUMBER OF FREE MASSES =  %d ------------\n", NMF);
    int_t NX = 6*NMF;
    int_t NU = 3;
    int_t jj;

    // Problem data
    int_t   N                   = NN;
    real_t  *x0;
    real_t  *w;  // States and controls stacked
    real_t  *Q;
    real_t  *R;
    real_t  *xref;
    real_t  *uref;
    int_t   max_sqp_iters       = 20;
    int_t   max_iters           = 1;
    real_t  *x_end;
    real_t  *u_end;
    FILE *initStates, *refStates;
    real_t feas, stepX, stepU;

    d_zeros(&x0, NX, 1);
    d_zeros(&xref, NX, 1);
    d_zeros(&uref, NU, 1);
    d_zeros(&x_end, NX, 1);
    d_zeros(&u_end, NU, 1);
    d_zeros(&w, NN*(NX+NU)+NX, 1);
    d_zeros(&Q, NX, NX);
    d_zeros(&R, NU, NU);

    switch (NMF) {
    case 1:
        initStates = fopen("../test/casadi_chain/x0_nm2.txt", "r");
        break;
    case 2:
        initStates = fopen("../test/casadi_chain/x0_nm3.txt", "r");
        break;
    case 3:
        initStates = fopen("../test/casadi_chain/x0_nm4.txt", "r");
        break;
    case 4:
        initStates = fopen("../test/casadi_chain/x0_nm5.txt", "r");
        break;
    case 5:
        initStates = fopen("../test/casadi_chain/x0_nm6.txt", "r");
        break;
    case 6:
        initStates = fopen("../test/casadi_chain/x0_nm7.txt", "r");
        break;
    case 7:
        initStates = fopen("../test/casadi_chain/x0_nm8.txt", "r");
        break;
    default:
        initStates = fopen("../test/casadi_chain/x0_nm9.txt", "r");
        break;
    }
    for (int_t i = 0; i < NX; i++) {
        nil = fscanf(initStates, "%lf", &x0[i]);
    }
    fclose(initStates);

    switch (NMF) {
    case 1:
        refStates = fopen("../test/casadi_chain/xN_nm2.txt", "r");
        break;
    case 2:
        refStates = fopen("../test/casadi_chain/xN_nm3.txt", "r");
        break;
    case 3:
        refStates = fopen("../test/casadi_chain/xN_nm4.txt", "r");
        break;
    case 4:
        refStates = fopen("../test/casadi_chain/xN_nm5.txt", "r");
        break;
    case 5:
        refStates = fopen("../test/casadi_chain/xN_nm6.txt", "r");
        break;
    case 6:
        refStates = fopen("../test/casadi_chain/xN_nm7.txt", "r");
        break;
    case 7:
        refStates = fopen("../test/casadi_chain/xN_nm8.txt", "r");
        break;
    default:
        refStates = fopen("../test/casadi_chain/xN_nm9.txt", "r");
        break;
    }
    for (int_t i = 0; i < NX; i++) {
        nil = fscanf(refStates, "%lf", &xref[i]);
    }
    fclose(refStates);

    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = xref[j];
    }
    for (int_t j = 0; j < NX; j++) w[NN*(NX+NU)+j] = xref[j];

    for (int_t i = 0; i < NX; i++) Q[i*(NX+1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i*(NU+1)] = 0.01;

    // Integrator structs
    real_t T = 0.2;
    sim_in  sim_in;
    sim_out sim_out;
    sim_in.nSteps = 10;
    sim_in.step = T/sim_in.nSteps;
    sim_in.nx = NX;
    sim_in.nu = NU;

    switch (NMF) {
    case 1:
        sim_in.VDE_fun = &VDE_fun_nm2;
        break;
    case 2:
        sim_in.VDE_fun = &VDE_fun_nm3;
        break;
    case 3:
        sim_in.VDE_fun = &VDE_fun_nm4;
        break;
    case 4:
        sim_in.VDE_fun = &VDE_fun_nm5;
        break;
    case 5:
        sim_in.VDE_fun = &VDE_fun_nm6;
        break;
    case 6:
        sim_in.VDE_fun = &VDE_fun_nm7;
        break;
    case 7:
        sim_in.VDE_fun = &VDE_fun_nm8;
        break;
    default:
        sim_in.VDE_fun = &VDE_fun_nm9;
        break;
    }

    sim_in.x = malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = malloc(sizeof(*sim_in.u) * (NU));

    sim_out.xn = malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.Sx = malloc(sizeof(*sim_out.Sx) * (NX*NX));
    sim_out.Su = malloc(sizeof(*sim_out.Su) * (NX*NU));

    sim_erk_workspace erk_work;
    sim_RK_opts rk_opts;
    sim_erk_create_opts(4, &rk_opts);
    sim_erk_create_workspace(&sim_in, &rk_opts, &erk_work);

    int_t nx[NN+1] = {0};
    int_t nu[NN] = {0};
    int_t nb[NN+1] = {0};
    int_t nc[NN+1] = {0};
    for (int_t i = 0; i < N; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[N] = NX;

    /************************************************
    * box constraints
    ************************************************/

    int *idxb0;
    i_zeros(&idxb0, NX+NU, 1);
    real_t *lb0;
    d_zeros(&lb0, NX+NU, 1);
    real_t *ub0;
    d_zeros(&ub0, NX+NU, 1);
    for (jj = 0; jj < NX; jj++) {
        lb0[jj] = x0[jj];   //   xmin
        ub0[jj] = x0[jj];   //   xmax
        idxb0[jj] = jj;
    }
    for (; jj < NX+NU; jj++) {
        lb0[jj] = -UMAX;  //   umin
        ub0[jj] = UMAX;   //   umax
        idxb0[jj] = jj;
    }
    nb[0] = NX+NU;

    int *idxb1;
    i_zeros(&idxb1, NU, 1);
    double *lb1[N-1];
    double *ub1[N-1];
    for (int_t i = 0; i < N-1; i++) {
    d_zeros(&lb1[i], NU, 1);
    d_zeros(&ub1[i], NU, 1);
//    for (jj = 0; jj < nbx; jj++) {
//        lb1[jj] = -4.0;  //   xmin
//        ub1[jj] = +4.0;   //   xmax
//        idxb1[jj] = jj;
//    }
    for (jj = 0; jj < NU; jj++) {
        lb1[i][jj] = -UMAX;  //   umin
        ub1[i][jj] = UMAX;   //   umax
        idxb1[jj] = NX+jj;
    }
    }

    real_t *pA[N];
    real_t *pB[N];
    real_t *pb[N];
    real_t *pQ[N+1];
    real_t *pS[N];
    real_t *pR[N];
    real_t *pq[N+1];
    real_t *pr[N];
    real_t *px[N+1];
    real_t *pu[N];
    for (int_t i = 0; i < N; i++) {
        d_zeros(&pA[i], nx[i+1], nx[i]);
        d_zeros(&pB[i], nx[i+1], nu[i]);
        d_zeros(&pb[i], nx[i+1], 1);
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
    }
    d_zeros(&pq[N], nx[N], 1);
    d_zeros(&px[N], nx[N], 1);

    real_t *hlb[N+1];
    real_t *hub[N+1];
    int *hidxb[N+1];

    hlb[0] = lb0;
    hub[0] = ub0;
    hidxb[0] = idxb0;
    for (int_t i = 1; i < N; i++) {
        hlb[i] = lb1[i-1];
        hub[i] = ub1[i-1];
        hidxb[i] = idxb1;
        nb[i] = NU;
    }

    // Allocate OCP QP variables
    ocp_qp_input qp_in;
    qp_in.N = N;
    ocp_qp_output qp_out;
    ocp_qp_condensing_qpoases_args args;
    real_t *work = NULL;
    qp_in.nx = nx;
    qp_in.nu = nu;
    qp_in.nb = nb;
    qp_in.nc = nc;
    for (int_t i = 0; i < N; i++) {
        pQ[i] = Q;
        pR[i] = R;
    }
    pQ[N] = Q;
    qp_in.Q = (const real_t **) pQ;
    qp_in.S = (const real_t **) pS;
    qp_in.R = (const real_t **) pR;
    qp_in.q = (const real_t **) pq;
    qp_in.r = (const real_t **) pr;
    qp_in.A = (const real_t **) pA;
    qp_in.B = (const real_t **) pB;
    qp_in.b = (const real_t **) pb;
    qp_in.lb = (const real_t **) hlb;
    qp_in.ub = (const real_t **) hub;
    qp_in.idxb = (const int **) hidxb;
    qp_out.x = px;
    qp_out.u = pu;

    initialise_qpoases(&qp_in);

    acado_timer timer;
    real_t timings = 0;
//    for (int_t iter = 0; iter < max_iters; iter++) {
//        printf("\n------ TIME STEP %d ------\n", iter);

        acado_tic(&timer);
        for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
            printf("--- ITERATION %d, ", sqp_iter);

            feas = -1e10; stepX = -1e10; stepU = -1e10;

            for (int_t i = 0; i < N; i++) {
                // Pass state and control to integrator
                for (int_t j = 0; j < NX; j++) sim_in.x[j] = w[i*(NX+NU)+j];
                for (int_t j = 0; j < NU; j++) sim_in.u[j] = w[i*(NX+NU)+NX+j];
                integrate(&sim_in, &sim_out, &rk_opts, &erk_work);

//                #ifdef DEBUG
//                print_matrix("stdout", sim_out.xn, 1, NX);
//                print_matrix("stdout", sim_out.Sx, NX, NX);
//                print_matrix("stdout", sim_out.Su, NX, NU);
//                #endif  // DEBUG

                // Update bounds:
                if ( i == 0 ) {
                    for (int_t j = 0; j < NU; j++) {
                        lb0[NX+j] = -UMAX - w[NX+j];
                        ub0[NX+j] = UMAX - w[NX+j];
                    }
                } else {
                    for (int_t j = 0; j < NU; j++) {
                        lb1[i-1][j] = -UMAX - w[i*(NX+NU)+NX+j];
                        ub1[i-1][j] = UMAX - w[i*(NX+NU)+NX+j];
                    }
                }

                // Construct QP matrices
                for (int_t j = 0; j < NX; j++) {
                    pq[i][j] = Q[j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
                }
                for (int_t j = 0; j < NU; j++) {
                    pr[i][j] = R[j*(NU+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
                }
                for (int_t j = 0; j < NX; j++) {
                    pb[i][j] = sim_out.xn[j] - w[(i+1)*(NX+NU)+j];
                    if (fabs(pb[i][j]) > feas) feas = fabs(pb[i][j]);
                    for (int_t k = 0; k < NX; k++)
                        pA[i][j*NX+k] = sim_out.Sx[j*NX+k];  // COLUMN MAJOR FROM CASADI
                }
                for (int_t j = 0; j < NU; j++)
                    for (int_t k = 0; k < NX; k++)
                        pB[i][j*NX+k] = sim_out.Su[j*NX+k];  // COLUMN MAJOR FROM CASADI
            }
            for (int_t j = 0; j < NX; j++) {
                lb0[j] = (x0[j]-w[j]);
            }
            for (int_t j = 0; j < NX; j++) {
                pq[N][j] = Q[j]*(w[N*(NX+NU)+j]-xref[j]);
            }

            // Set updated bounds:
            hlb[0] = lb0;
            hub[0] = ub0;
            for (int_t i = 1; i < N; i++) {
                hlb[i] = lb1[i-1];
                hub[i] = ub1[i-1];
            }

            int status = 0;
            status = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, work);
            if (status) {
                printf("qpOASES returned error status %d\n", status);
                return -1;
            }
            for (int_t i = 0; i < N; i++) {
                for (int_t j = 0; j < NX; j++) {
                    w[i*(NX+NU)+j] += qp_out.x[i][j];
                    if (fabs(qp_out.x[i][j]) > stepX) stepX = fabs(qp_out.x[i][j]);
                }
                for (int_t j = 0; j < NU; j++) {
                    w[i*(NX+NU)+NX+j] += qp_out.u[i][j];
                    if (fabs(qp_out.u[i][j]) > stepU) stepU = fabs(qp_out.u[i][j]);
                }
            }
            for (int_t j = 0; j < NX; j++) {
                w[N*(NX+NU)+j] += qp_out.x[N][j];
                if (fabs(qp_out.x[N][j]) > stepX) stepX = fabs(qp_out.x[N][j]);
            }

            fprintf(stdout, " Infeasibility: %+.3e , step X: %+.3e, step U: %+.3e \n",
                    feas, stepX, stepU);
        }
//        for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
//        shift_states(w, x_end, N);
//        shift_controls(w, u_end, N);
        timings += acado_toc(&timer);
//    }

    printf("\nAverage of %.3f ms per iteration.\n\n", 1e3*timings/(max_sqp_iters*max_iters));

    #ifdef DEBUG
    print_matrix("stdout", w, NX+NU, N);
    #endif  // DEBUG
    }
    return 0*nil;
}

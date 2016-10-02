#include <stdio.h>
#include <stdlib.h>
#if defined(__APPLE__)
#include <mach/mach_time.h>
#else
#include <sys/stat.h>
#include <sys/time.h>
#endif
#include "acados/utils/types.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "external/hpmpc/include/aux_d.h"
#include "acados/model.h"

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

#ifdef DEBUG
static void print_states_controls(real_t *w) {
    int_t   N = NNN;
    printf("node\tx\t\t\t\tu\n");
    for (int_t i = 0; i < N; i++) {
        printf("%4d\t%+e %+e\t%+e\n", i, w[i*(NX+NU)], w[i*(NX+NU)+1], w[i*(NX+NU)+2]);
    }
    printf("%4d\t%+e %+e\n", N, w[N*(NX+NU)], w[N*(NX+NU)+1]);
}
#endif  // DEBUG

static void shift_states(real_t *w, real_t *x_end) {
    for (int_t i = 0; i < NNN; i++) {
        for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = w[(i+1)*(NX+NU)+j];
    }
    for (int_t j = 0; j < NX; j++) w[NNN*(NX+NU)+j] = x_end[j];
}

static void shift_controls(real_t *w, real_t *u_end) {
    for (int_t i = 0; i < NNN-1; i++) {
        for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = w[(i+1)*(NX+NU)+NX+j];
    }
    for (int_t j = 0; j < NU; j++) w[(NNN-1)*(NX+NU)+NX+j] = u_end[j];
}

// Simple SQP example for acados
int main() {
    // Problem data
    int_t   N               = NNN;
    real_t  x0[NX]          = {0.5, 0};      // Initial states
    real_t  w[NNN*(NX+NU)+NX] = {0};        // States and controls stacked
    real_t  Q[NX*NX]        = {0};          // Weight on the states
    real_t  R[NU*NU]        = {0};          // Weight on the controls
    real_t  xref[NX]        = {0};          // State reference
    real_t  uref[NX]        = {0};          // Control reference
    int_t   max_sqp_iters   = 1;           // Number of SQP iterations
    int_t   max_iters       = 10000;      // Number of NMPC samples
    real_t  x_end[NX]       = {0};
    real_t  u_end[NU]       = {0};

    for (int_t i = 0; i < NX; i++) Q[i*(NX+1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i*(NU+1)] = 0.05;

    // Integrator structs
    real_t T = 0.1;
    sim_in  sim_in;
    sim_out sim_out;
    sim_in.nSteps = 10;
    sim_in.step = T/sim_in.nSteps;
    sim_in.VDE_fun = &VDE_fun;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.x = malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = malloc(sizeof(*sim_in.u) * (NU));

    sim_out.xn = malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.Sx = malloc(sizeof(*sim_out.Sx) * (NX*NX));
    sim_out.Su = malloc(sizeof(*sim_out.Su) * (NX*NU));

    sim_erk_workspace erk_work;
    sim_RK_opts rk_opts;
    sim_erk_create_opts(4, &rk_opts);
    sim_erk_create_workspace(&sim_in, &rk_opts, &erk_work);

    int_t nx[NNN+1] = {0};
    int_t nu[NNN] = {0};
    int_t nb[NNN+1] = {0};
    int_t nc[NNN+1] = {0};
    for (int_t i = 0; i < N; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[N] = NX;

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
    real_t *px0[1];
    d_zeros(&px0[0], nx[0], 1);
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
    qp_in.lb = (const real_t **) px0;
    qp_out.x = px;
    qp_out.u = pu;

    initialise_qpoases(&qp_in);

    acado_timer timer;
    real_t timings = 0;
    for (int_t iter = 0; iter < max_iters; iter++) {
        // printf("\n------ ITERATION %d ------\n", iter);
        acado_tic(&timer);
        for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
            for (int_t i = 0; i < N; i++) {
                // Pass state and control to integrator
                for (int_t j = 0; j < NX; j++) sim_in.x[j] = w[i*(NX+NU)+j];
                for (int_t j = 0; j < NU; j++) sim_in.u[j] = w[i*(NX+NU)+NX+j];
                integrate(&sim_in, &sim_out, &rk_opts, &erk_work);
                // Construct QP matrices
                for (int_t j = 0; j < NX; j++) {
                    pq[i][j] = Q[j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
                }
                for (int_t j = 0; j < NU; j++) {
                    pr[i][j] = R[j*(NX+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
                }
                for (int_t j = 0; j < NX; j++) {
                    pb[i][j] = sim_out.xn[j] - w[(i+1)*(NX+NU)+j];
                    for (int_t k = 0; k < NX; k++) pA[i][j*NX+k] = sim_out.Sx[k*NX+j];
                }
                for (int_t j = 0; j < NU; j++)
                    for (int_t k = 0; k < NX; k++) pB[i][j*NX+k] = sim_out.Su[k*NU+j];
            }
            for (int_t j = 0; j < NX; j++) {
                px0[0][j] = (x0[j]-w[j]);
            }
            for (int_t j = 0; j < NX; j++) {
                pq[N][j] = Q[j]*(w[N*(NX+NU)+j]-xref[j]);
            }
            int status = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, work);
            if (status) {
                printf("qpOASES returned error status %d\n", status);
                return -1;
            }
            for (int_t i = 0; i < N; i++) {
                for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] += qp_out.x[i][j];
                for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] += qp_out.u[i][j];
            }
            for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] += qp_out.x[N][j];
        }
        for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
        shift_states(w, x_end);
        shift_controls(w, u_end);
        timings += acado_toc(&timer);
    }
    #ifdef DEBUG
    print_states_controls(&w[0]);
    #endif  // DEBUG
    printf("Average of %.3f ms per iteration.\n", 1e3*timings/max_iters);

    return 0;
}

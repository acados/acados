#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "acados/acados_types.h"
#include "acados/erk_integrator.h"
#include "acados/ocp_qp_condensing_qpoases.h"
#include "hpmpc/include/aux_d.h"

#include "unistd.h"
#include <mach/mach_time.h>

/** A structure for keeping internal timer data. */
typedef struct acado_timer_
{
	uint64_t tic;
	uint64_t toc;
	mach_timebase_info_data_t tinfo;
} acado_timer;

static void acado_tic( acado_timer* t )
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

static real_t acado_toc( acado_timer* t )
{

    uint64_t duration; /* elapsed time in clock cycles*/

    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (real_t)duration / 1e9;
}

// static void print_states_controls(real_t *w) {
//     int_t   N = NNN;
//     printf("node\tx\t\t\t\tu\n");
//     for (int_t i = 0; i < N; i++) {
//         printf("%4d\t%+e %+e\t%+e\n", i, w[i*(NX+NU)],w[i*(NX+NU)+1],w[i*(NX+NU)+2]);
//     }
//     printf("%4d\t%+e %+e\n", N, w[N*(NX+NU)],w[N*(NX+NU)+1]);
// }

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
    real_t  x0[NX]          = {0.5,0};      // Initial states
    real_t  w[NNN*(NX+NU)+NX] = {0};        // States and controls stacked
    real_t  Q[NX*NX]        = {0};          // Weight on the states
    real_t  R[NU*NU]        = {0};          // Weight on the controls
    real_t  xref[NX]        = {0};          // State reference
    real_t  uref[NX]        = {0};          // Control reference
    int_t   max_sqp_iters   = 1;           // Number of SQP iterations
    int_t   max_iters       = 10;           // Number of NMPC samples
    real_t  x_end[NX]       = {0};
    real_t  u_end[NU]       = {0};

    for (int_t i = 0; i < NX; i++) Q[i*(NX+1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i*(NU+1)] = 0.05;

    // Integrator structs
    real_t T = 0.1;
    sim_in  in;
    sim_out out;
    in.nSteps = 10;
    in.step = T/in.nSteps;

    int_t nx[NNN+1] = {0};
    int_t nu[NNN] = {0};
    int_t nb[NNN+1] = {0};
    int_t nc[NNN+1] = {0};
    for (int_t i = 0; i < N; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[N] = NX;

    // Allocate OCP QP variables
    ocp_qp_input qp_in;
    qp_in.N = N;
    ocp_qp_output qp_out;
    ocp_qp_condensing_qpoases_args args;
    real_t *work = NULL;

    // Fill in dimensions
    qp_in.nx = nx;
    qp_in.nu = nu;
    qp_in.nb = nb;
    qp_in.nc = nc;

    qp_in.A = malloc(sizeof(*qp_in.A) * N);
    qp_in.B = malloc(sizeof(*qp_in.B) * N);
    qp_in.b = malloc(sizeof(*qp_in.b) * N);
    qp_in.Q = malloc(sizeof(*qp_in.Q) * (N+1));
    qp_in.S = malloc(sizeof(*qp_in.S) * N);
    qp_in.R = malloc(sizeof(*qp_in.R) * N);
    qp_in.q = malloc(sizeof(*qp_in.q) * (N+1));
    qp_in.r = malloc(sizeof(*qp_in.r) * N);
    qp_out.x = malloc(sizeof(*qp_out.x) * (N+1));
    qp_out.u = malloc(sizeof(*qp_out.u) * N);
    for (int_t i = 0; i < N; i++) {
        d_zeros(&qp_in.A[i], nx[i+1], nx[i]);
        d_zeros(&qp_in.B[i], nx[i+1], nu[i]);
        d_zeros(&qp_in.b[i], nx[i+1], 1);
        d_zeros(&qp_in.S[i], nu[i], nx[i]);
        d_zeros(&qp_in.q[i], nx[i], 1);
        d_zeros(&qp_in.r[i], nu[i], 1);
        d_zeros(&qp_out.x[i], nx[i], 1);
        d_zeros(&qp_out.u[i], nu[i], 1);
    }
    d_zeros(&qp_in.q[N], nx[N], 1);
    d_zeros(&qp_out.x[N], nx[N], 1);

    // Fill in objective
    for (int_t i = 0; i < N; i++) {
        qp_in.Q[i] = Q;
        qp_in.R[i] = R;
    }
    qp_in.Q[N] = Q;

    qp_in.lb = malloc(sizeof(*qp_in.lb));
    d_zeros(&qp_in.lb[0], nx[0], 1);

    initialise_qpoases(&qp_in);

    real_t timings[max_iters];
    acado_timer timer;
    for (int_t iter = 0; iter < max_iters; iter++) {
        // printf("\n------ ITERATION %d ------\n", iter);
        acado_tic(&timer);
        for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
            for (int_t i = 0; i < N; i++) {
                // Pass state and control to integrator
                for (int_t j = 0; j < NX; j++) in.x[j] = w[i*(NX+NU)+j];
                for (int_t j = 0; j < NU; j++) in.u[j] = w[i*(NX+NU)+NX+j];
                integrate(&in, &out);
                // Construct QP matrices
                for (int_t j = 0; j < NX; j++) {
                    qp_in.q[i][j] = Q[j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
                }
                for (int_t j = 0; j < NU; j++) {
                    qp_in.r[i][j] = R[j*(NX+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
                }
                for (int_t j = 0; j < NX; j++) {
                    qp_in.b[i][j] = out.xn[j] - w[(i+1)*(NX+NU)+j];
                    for (int_t k = 0; k < NX; k++) qp_in.A[i][j*NX+k] = out.Sx[k*NX+j];
                }
                for (int_t j = 0; j < NU; j++)
                    for (int_t k = 0; k < NX; k++) qp_in.B[i][j*NX+k] = out.Su[k*NU+j];
            }
            for (int_t j = 0; j < NX; j++) {
                qp_in.lb[0][j] = (x0[j]-w[j]);
            }
            for (int_t j = 0; j < NX; j++) {
                qp_in.q[N][j] = Q[j]*(w[N*(NX+NU)+j]-xref[j]);
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
        // print_states_controls(&w[0]);
        for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
        shift_states(w, x_end);
        shift_controls(w, u_end);
        timings[iter] = acado_toc(&timer);
    }
    for (int_t iter = 0; iter < max_iters; iter++) {
        printf("%d:\t%f\tms\n", iter, 1e3*timings[iter]);
    }

    return 0;
}

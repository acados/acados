/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// #define PLOT_RESULTS
// #define FP_EXCEPTIONS

#ifdef PLOT_RESULTS
#define _GNU_SOURCE
#endif  // PLOT_RESULTS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// flush denormals to zero
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) ||  \
    defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || \
    defined(TARGET_AMD_SSE3)
#include <xmmintrin.h>  // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

#include "hpmpc/include/aux_d.h"
#include "hpmpc/include/lqcp_solvers.h"

#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/utils/timing.h"
#include "acados/utils/math.h"
#include "acados/utils/types.h"

// define IP solver arguments && number of repetitions
#define NREP 10
#define MAX_IP_ITER 30
#define TOL 1e-6
#define MINSTEP 1e-8

#define NN 100
#define MM 10
#define NX 4
#define NU 1
#define NBU 1
#define NBX 0  // TODO(Andrea): adding bounds gives MIN_STEP
#define NSIM 1
#define GAMMA 1e-2
#define UMAX 1

int vdeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

#ifdef DEBUG
static void print_states_controls(real_t *w, int_t N) {
    printf("node\tx\t\t\t\t\t\tu\n");
    for (int_t i = 0; i < N; i++) {
        printf("%4d\t%+e %+e %+e %+e\t%+e\n", i, w[i*(NX+NU)], w[i*(NX+NU)+1],
        w[i*(NX+NU)+2], w[i*(NX+NU)+3], w[i*(NX+NU)+4]);
    }
    printf("%4d\t%+e %+e %+e %+e \n", N, w[N*(NX+NU)],
    w[N*(NX+NU)+1], w[N*(NX+NU)+2], w[N*(NX+NU)+3]);
}
#endif  // DEBUG

#ifdef PLOT_RESULTS
// static void plot_states_controls(real_t *w, int_t nx, int_t nu, int_t N, real_t T ) {
static void plot_states_controls(real_t *w, real_t T) {
      double t_grid[NN];
      for (int_t i = 0; i < NN; i++) t_grid[i] = i*T;

      FILE *gnuplotPipe, *tempDataFile;
      char *x1_temp_file;
      char *x2_temp_file;
      char *u1_temp_file;
      double x, y;
      int i;
      x1_temp_file = "x1";
      gnuplotPipe = popen("gnuplot -persist", "w");
      if (gnuplotPipe) {
          fprintf(gnuplotPipe, "set multiplot layout 3,1\n");

          // Plot x1
          tempDataFile = fopen(x1_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x1_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot x2
          x2_temp_file = "x2";
          tempDataFile = fopen(x2_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x2_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+1];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u1
          u1_temp_file = "u1";
          tempDataFile = fopen(u1_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with steps lt rgb \"red\" \n", u1_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+4];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          printf("Press enter to continue...");
          getchar();
          remove(x1_temp_file);
          remove(x2_temp_file);
          remove(u1_temp_file);

          fprintf(gnuplotPipe, "exit gnuplot\n");
      } else {
          printf("gnuplot not found...");
      }
}

// static void plot_states_controls_cl(real_t *w, real_t T) {
//       double t_grid[NN];
//       for (int_t i = 0; i < NSIM; i++) t_grid[i] = i*T;
//
//       FILE *gnuplotPipe, *tempDataFile;
//       char *x1_temp_file;
//       char *x2_temp_file;
//       char *u1_temp_file;
//       double x, y;
//       int i;
//       x1_temp_file = "x1";
//       gnuplotPipe = popen("gnuplot -persist", "w");
//       if (gnuplotPipe) {
//           fprintf(gnuplotPipe, "set multiplot layout 3,1\n");
//
//           // Plot x1
//           tempDataFile = fopen(x1_temp_file, "w");
//           fprintf(gnuplotPipe, "set grid ytics\n");
//           fprintf(gnuplotPipe, "set grid xtics\n");
//           fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
//           fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x1_temp_file);
//           fflush(gnuplotPipe);
//           for (i=0; i < NSIM; i++) {
//               x = t_grid[i];
//               y = w[i*(NX+NU)];
//               fprintf(tempDataFile, "%lf %lf\n", x, y);
//           }
//           fclose(tempDataFile);
//
//           // Plot x2
//           x2_temp_file = "x2";
//           tempDataFile = fopen(x2_temp_file, "w");
//           fprintf(gnuplotPipe, "set grid ytics\n");
//           fprintf(gnuplotPipe, "set grid xtics\n");
//           fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
//           fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x2_temp_file);
//           fflush(gnuplotPipe);
//           for (i=0; i < NSIM; i++) {
//               x = t_grid[i];
//               y = w[i*(NX+NU)+1];
//               fprintf(tempDataFile, "%lf %lf\n", x, y);
//           }
//           fclose(tempDataFile);
//
//           // Plot u1
//           u1_temp_file = "u1";
//           tempDataFile = fopen(u1_temp_file, "w");
//           fprintf(gnuplotPipe, "set grid ytics\n");
//           fprintf(gnuplotPipe, "set grid xtics\n");
//           fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
//           fprintf(gnuplotPipe, "plot \"%s\" with steps lt rgb \"red\" \n", u1_temp_file);
//           fflush(gnuplotPipe);
//           for (i=0; i < NSIM; i++) {
//               x = t_grid[i];
//               y = w[i*(NX+NU)+4];
//               fprintf(tempDataFile, "%lf %lf\n", x, y);
//           }
//           fclose(tempDataFile);
//
//           printf("Press enter to continue...");
//           getchar();
//           remove(x1_temp_file);
//           remove(x2_temp_file);
//           remove(u1_temp_file);
//
//           fprintf(gnuplotPipe, "exit gnuplot\n");
//       } else {
//           printf("gnuplot not found...");
//       }
// }
#endif  // PLOT_RESULTS


// static void shift_states(real_t *w, real_t *x_end, int_t N) {
//     for (int_t i = 0; i < N; i++) {
//         for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = w[(i+1)*(NX+NU)+j];
//     }
//     for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] = x_end[j];
// }

// static void shift_controls(real_t *w, real_t *u_end, int_t N) {
//     for (int_t i = 0; i < N-1; i++) {
//         for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = w[(i+1)*(NX+NU)+NX+j];
//     }
//     for (int_t j = 0; j < NU; j++) w[(N-1)*(NX+NU)+NX+j] = u_end[j];
// }

// TODO(Andrea): need to fix this stuff below...

extern int d_ip2_res_mpc_hard_work_space_size_bytes_libstr(int N, int *nx,
  int *nu, int *nb, int *ng);

extern int d_back_ric_rec_work_space_size_bytes_libstr(int N, int *nx, int *nu,
    int *nb, int *ng);  // TODO(Andrea): ask Gian how to fix this

extern  int d_size_strmat(int m, int n);
extern  int d_size_strvec(int m);
// extern FILE *popen(char *command, const char *type);

// Simple SQP example for acados
int main() {
    #ifdef FP_EXCEPTIONS
      feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
      // feenableexcept(FE_INVALID);
    #endif  // FP_EXCEPTIONS

    // Problem data
    int_t   N                         = NN;
    int_t   M                         = MM;
    real_t  x0[NX]                    = {0.0, 0.1, 0.0, 0.0};
    real_t  w[NN*(NX+NU)+NX]          = {0};  // States and controls stacked
    // real_t  w_cl[NSIM*(NX+NU)]        = {0};  // States and controls stacked closed loop
    // real_t  pi_n[NN*(NX)]             = {0};
    real_t  t_n[2*((NBX+NBU)*NN + NBX)]   = {0};
    real_t  lam_n[2*((NBX+NBU)*NN + NBX)] = {0};
    real_t  Q[NX*NX]                  = {0};
    real_t  R[NU*NU]                  = {0};
    real_t  xref[NX]                  = {0};
    real_t  uref[NX]                  = {0};
    real_t  lam_init                  = {0.0001};
    real_t  t_init                    = {1};
    real_t sigma_mu                   = {0.001};
    // real_t  pi_init                   = {0.1};
    // int_t   qp_iters               = 1;
    // int_t   max_iters               = 100;
    // real_t  x_min[NBX]              = {-10, -10, -10, -10};
    real_t  x_min[NBX]                = {};
    // real_t  x_max[NBX]             = {10, 10, 10, 10};
    real_t  x_max[NBX]                = {};
    real_t  u_min[NBU]                = {-UMAX};
    real_t  u_max[NBU]                = {UMAX};

    // Q[0*(NX+1)] = 1e-3;
    // Q[1*(NX+1)] = 5e-1;
    // Q[2*(NX+1)] = 1.0;
    // Q[3*(NX+1)] = 2e-3;

    Q[0*(NX+1)] = 1e-3;
    Q[1*(NX+1)] = 5e1;
    Q[2*(NX+1)] = 0.01;
    Q[3*(NX+1)] = 2e-3;

    R[0*(NU+1)] = 1e-4;

    // double sigma_mu = 0.01;

    // Integrator structs
    real_t T = 0.01;
    sim_in  sim_in;
    sim_out sim_out;
    sim_in.num_steps = 1;
    sim_in.step = T/sim_in.num_steps;
    sim_in.vde = &vdeFun;
    sim_in.forward_vde_wrapper = &vde_fun;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.sens_forw = true;
    sim_in.sens_adj = false;
    sim_in.sens_hess = false;
    sim_in.num_forw_sens = NX+NU;

    sim_in.x = malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = malloc(sizeof(*sim_in.u) * (NU));
    sim_in.S_forw = malloc(sizeof(*sim_in.S_forw) * (NX*(NX+NU)));
    for (int_t i = 0; i < NX*(NX+NU); i++) sim_in.S_forw[i] = 0.0;
    for (int_t i = 0; i < NX; i++) sim_in.S_forw[i*(NX+1)] = 1.0;

    sim_out.xn = malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.S_forw = malloc(sizeof(*sim_out.S_forw) * (NX*(NX+NU)));

    sim_info erk_info;
    sim_out.info = &erk_info;

    void *erk_work;
    sim_rk_opts rk_opts;
    sim_erk_create_arguments(&rk_opts, 4);
    int_t sim_workspace_size = sim_erk_calculate_workspace_size(&sim_in, &rk_opts);
    erk_work = (void *) malloc(sim_workspace_size);
    // sim_erk_create_workspace(&sim_in, &rk_opts, &erk_work);
// =======
//     sim_rk_opts rk_opts;
//     sim_erk_create_arguments(&rk_opts, 4);
//     void *erk_work;
//     int_t erk_workspace_size = sim_erk_calculate_workspace_size(&sim_in, &rk_opts);
//     erk_work = (void *) malloc(erk_workspace_size);
// >>>>>>> 33c021b0a2ec9d74225c1c5ea688af5541ca3319

    int_t nx[NN+1] = {0};
    int_t nu[NN+1] = {0};
    int_t nb[NN+1] = {0};
    int_t nc[NN+1] = {0};
    for (int_t i = 0; i < N; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[0] = 0;  // x0 is eliminated
    nx[N] = NX;
    nu[N] = 0;

    real_t *pA[N];
    real_t *pB[N];
    real_t *pb[N];
    real_t *pQ[N+1];
    real_t *pS[N+1];
    real_t *pR[N+1];
    real_t *pq[N+1];
    real_t *pr[N+1];
    real_t *px[N+1];
    real_t *pu[N+1];
    real_t *px0[1];
    int    *hidxb[N+1];
    real_t *pub[N+1];
    real_t *plb[N+1];
    real_t *pC[N+1];
    real_t *pD[N+1];
    real_t *plg[N+1];
    real_t *pug[N+1];

    /************************************************
    * box constraints
    ************************************************/
    int ii, jj;
    nb[0] = NBU;
    for (ii = 1; ii < N; ii++ ) nb[ii] = NBU + NBX;
    nb[N] = NBX;

    // int *idxb0;
    // int_zeros(&idxb0, nb[0], 1);
    // for (jj = 0; jj < NBX+NBU; jj++ ) idxb0[jj] = jj;
    //
    // int *idxb1;
    // int_zeros(&idxb1, nb[1], 1);
    // for (jj = 0; jj < NBX; jj++ ) idxb1[jj] = jj;
    // for (; jj < NBX+NBU; jj++ ) idxb1[jj] = NX+jj;
    //
    // int *idxbN;
    // int_zeros(&idxbN, nb[N], 1);
    // for (jj = 0; jj < NBX; jj++ ) idxbN[jj] = jj;
    // for (; jj < NBX+NBU; jj++ ) idxbN[jj] = NX+jj;
    //
    // hidxb[0] = idxb0;
    // for (ii = 1; ii < N; ii++ ) hidxb[ii] = idxb1;
    // hidxb[N] = idxbN;

    int *idxb0;  // TODO(Andrea): need to swap these guys in the interface...
    int_zeros(&idxb0, nb[0], 1);
    for (jj = 0; jj < NBU; jj++ ) idxb0[jj] = jj;

    int *idxb1;
    int_zeros(&idxb1, nb[1], 1);
    for (jj = 0; jj < NBU; jj++ ) idxb1[jj] = jj;
    for (; jj < NBX+NBU; jj++ ) idxb1[jj] = NU+jj;

    int *idxbN;
    int_zeros(&idxbN, nb[N], 1);
    // for (jj = 0; jj < NBU; jj++ ) idxbN[jj] = jj;
    for ( jj = 0; jj < NBX; jj++ ) idxbN[jj] = jj;

    hidxb[0] = idxb0;
    for (ii = 1; ii < N; ii++ ) hidxb[ii] = idxb1;
    // for (ii = M; ii < N; ii++ ) hidxb[ii] = idxb_tight;
    hidxb[N] = idxbN;

    d_zeros(&px0[0], nx[0], 1);
    d_zeros(&plb[0], NBU, 1);
    d_zeros(&pub[0], NBU, 1);
    for (int_t i = 0; i < N; i++) {
        d_zeros(&pA[i], nx[i+1], nx[i]);
        d_zeros(&pB[i], nx[i+1], nu[i]);
        d_zeros(&pb[i], nx[i+1], 1);
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
        if (i >= 1) d_zeros(&plb[i], NBU+NBX, 1);
        if (i >= 1) d_zeros(&pub[i], NBU+NBX, 1);
    }
    d_zeros(&plb[N], NBX, 1);
    d_zeros(&pub[N], NBX, 1);
    d_zeros(&pq[N], nx[N], 1);
    d_zeros(&px[N], nx[N], 1);
    // hidxb[N] = idxbN;

    //    d_print_mat(nx, 1, b, nx);
    //    d_print_mat(nx, 1, b0, nx);

    // then A0 is a matrix of size 0x0
    double *A0;
    d_zeros(&A0, 0, 0);

    int ng = 0;  // 4;  // number of general constraints
    int ngN = 0;  // 4;  // number of general constraints at the last stage

    int ngg[N + 1];
    for (ii = 0; ii < N; ii++) ngg[ii] = ng;
    ngg[N] = ngN;

    /************************************************
    * general constraints
    ************************************************/

    double *C0;
    d_zeros(&C0, 0, NX);
    double *D0;
    d_zeros(&D0, 0, NU);
    double *lg0;
    d_zeros(&lg0, 0, 1);
    double *ug0;
    d_zeros(&ug0, 0, 1);

    double *C;
    d_zeros(&C, 0, NX);
    double *D;
    d_zeros(&D, 0, NU);
    double *lg;
    d_zeros(&lg, 0, 1);
    double *ug;
    d_zeros(&ug, 0, 1);

    double *CN;
    d_zeros(&CN, ngN, NX);
    double *DN;
    d_zeros(&DN, ngN, nu[N]);
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

    pC[0] = C0;
    pD[0] = D0;
    plg[0] = lg0;
    pug[0] = ug0;
    real_t *ppi[N];
    real_t *plam[N+1];

    double *pt[N+1];

    double *lam_in[N+1];
    double *t_in[N+1];
    double *ux_in[N+1];

    ii = 0;
    d_zeros(&ppi[ii], nx[ii+1], 1);
    d_zeros(&plam[ii], 2*nb[ii]+2*nb[ii], 1);
    for (ii = 0; ii < N; ii++) {
        pC[ii] = C;
        pD[ii] = D;
        plg[ii] = lg;
        pug[ii] = ug;

        d_zeros(&ppi[ii], nx[ii+1], 1);
        d_zeros(&plam[ii], 2*nb[ii]+2*nb[ii], 1);
        d_zeros(&pt[ii], 2*nb[ii]+2*ngg[ii], 1);

        d_zeros(&lam_in[ii], 2*nb[ii]+2*ngg[ii], 1);
        d_zeros(&t_in[ii], 2*nb[ii]+2*ngg[ii], 1);
        d_zeros(&ux_in[ii], nx[ii]+nu[ii], 1);
    }

    d_zeros(&ppi[N], nx[N], 1);
    d_zeros(&plam[N], 2*nb[N]+2*nb[N], 1);
    d_zeros(&pt[N], 2*nb[N]+2*ngg[N], 1);

    d_zeros(&lam_in[N], 2*nb[N]+2*ngg[N], 1);
    d_zeros(&t_in[N], 2*nb[N]+2*ngg[N], 1);
    d_zeros(&ux_in[N], nx[N]+nu[N], 1);

    // // Init multipliers and slacks
    // for (jj = 0; jj < 2*nb[ii]+2*ngg[ii]; jj++) {
    //   lam_in[N][jj] = 1.0;
    //   t_in[N][jj] = 1.0;
    // }


    d_zeros(&plam[N], 2*nb[N]+2*nb[N], 1);

    pC[N] = CN;
    pD[N] = DN;
    plg[N] = lgN;
    pug[N] = ugN;

    /************************************************
    * solver arguments
    ************************************************/

    // solver arguments
    ocp_qp_hpmpc_args hpmpc_args;
    hpmpc_args.tol = TOL;
    hpmpc_args.max_iter = MAX_IP_ITER;
//  hpmpc_args.min_step = MINSTEP;
    hpmpc_args.mu0 = 1;
//  hpmpc_args.sigma_min = 1e-3;
    hpmpc_args.warm_start = 0;
    hpmpc_args.N2 = N;
    hpmpc_args.lam0 = lam_in;
    hpmpc_args.t0 = t_in;
    hpmpc_args.ux0 = ux_in;
    hpmpc_args.M = MM;
    hpmpc_args.sigma_mu = sigma_mu;

    // Allocate OCP QP variables
    ocp_qp_in qp_in;
    qp_in.N = N;
    ocp_qp_out qp_out;
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
    qp_in.lb = (const real_t **) plb;
    qp_in.ub = (const real_t **) pub;
    qp_in.idxb = (const int_t **) hidxb;
    qp_in.Cx = (const real_t **) pC;
    qp_in.Cu = (const real_t **) pD;
    qp_in.lc = (const real_t **) plg;
    qp_in.uc = (const real_t **) pug;

    qp_out.x = px;
    qp_out.u = pu;
    qp_out.pi = ppi;
    qp_out.lam = plam;
    qp_out.t = pt;

    void *workspace = 0;
    void *mem = 0;

    int_t work_space_size = 0;
    work_space_size = ocp_qp_hpmpc_calculate_workspace_size(&qp_in, &hpmpc_args);
    // printf("work_space_size = %i", work_space_size);
    v_zeros_align(&workspace, work_space_size);
    mem = ocp_qp_hpmpc_create_memory(&qp_in, &hpmpc_args);

    acados_timer timer;
    real_t timings = 0;
    real_t sum_timings = 0;
    real_t min_timings = 1000000;
    int_t status;

    // initialize nlp primal variables
    for (int_t i = 0; i < N; i++) {
      for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = 0.0;
      for (int_t j = 0; j < 2*(NBX+NBU); j++) lam_n[i*2*(NBX+NBU)+j]  = lam_init;
      for (int_t j = 0; j < 2*(NBX+NBU); j++) t_n[i*2*(NBX + NBU)+j]  = t_init;
      // for (int_t j = 0; j < NX; j++) pi_n[i*NX +j]  = pi_init;
      for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = 0.0;
    }

    for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] = 0.0;
    for (int_t j = 0; j < 2*NBX; j++) lam_n[N*2*(NBX+NBU)+j]  = lam_init;
    for (int_t j = 0; j < 2*NBX; j++) t_n[N*2*(NBX + NBU)+j]  = t_init;

    // initialize qp primal variables
    for (int_t j = 0; j < NU; j++) ux_in[0][j] = w[NX+j];
    for (int_t i = 1; i < N; i++) {
      for (int_t j = 0; j < NX; j++)  ux_in[i][j] = w[i*(NX+NU)+j];
      for (int_t j = 0; j < NU; j++)  ux_in[i][j] = w[i*(NX+NU)+NX+j];
    }

    for (int_t j = 0; j < NX; j++) ux_in[N][j] = w[N*(NX+NU)+j];

    for ( int_t ii = 0; ii < NX; ii++ ) w[ii] = x0[ii];

    for (int_t sim_iter = 0; sim_iter < NSIM; sim_iter++) {
      for (int_t iter = 0; iter < NREP; iter++) {
        acados_tic(&timer);
        // initialize nlp dual variables
        for (int_t i = M; i < N; i++) {
          for (int_t j  = 0; j < 2*nb[i]+2*ngg[i]; j++) {
            lam_in[i][j] = lam_n[2*(NBX + NBU)*i + j];
            t_in[i][j] = t_n[2*(NBX + NBU)*i + j];
          }
        }

        for (int_t j  = 0; j < 2*nb[N]+2*ngg[N]; j++) {
          lam_in[N][j] = lam_n[2*NBX*N + j];
          t_in[N][j] = t_n[2*NBX*N + j];
        }

        for (int_t i = 0; i < N; i++) {
            // Pass state and control to integrator
            for (int_t j = 0; j < NX; j++) sim_in.x[j] = w[i*(NX+NU)+j];
            for (int_t j = 0; j < NU; j++) sim_in.u[j] = w[i*(NX+NU)+NX+j];
            sim_erk(&sim_in, &sim_out, &rk_opts, 0, erk_work);
            // Construct QP matrices
            for (int_t j = 0; j < nx[i]; j++) {
                pq[i][j] = Q[j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
            }
            for (int_t j = 0; j < NU; j++) {
                pr[i][j] = R[j*(NU+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
            }
            for (int_t j = 0; j < NX; j++) {
                pb[i][j] = sim_out.xn[j] - w[(i+1)*(NX+NU)+j];
                for (int_t k = 0; k < nx[i]; k++) pA[i][j*NX+k] = sim_out.S_forw[j*(NX)+k];
            }
            for (int_t j = 0; j < NU; j++)
                for (int_t k = 0; k < NX; k++) pB[i][j*NX+k] = sim_out.S_forw[NX*NX + NX*j+k];


            for ( int_t j = 0; j < NBX; j++ ) plb[i][j+NBU] = x_min[j] - w[i*(NX+NU)+j];
            for ( int_t j = 0; j < NBX; j++ ) pub[i][j+NBU] = x_max[j] - w[i*(NX+NU)+j];
            for ( int_t j = 0; j < NBU; j++ ) plb[i][j] = u_min[j] - w[i*(NX+NU)+NX+j];
            for ( int_t j = 0; j < NBU; j++ ) pub[i][j] = u_max[j] - w[i*(NX+NU)+NX+j];
        }

        for ( int_t j = 0; j < NBX; j++ ) plb[N][j+NBU] = x_min[j] - w[N*(NX+NU)+j];
        for ( int_t j = 0; j < NBX; j++ ) pub[N][j+NBU] = x_max[j] - w[N*(NX+NU)+j];

        // for ( int_t j = 0; j < NBU; j++ ) plb[0][j] = u_min[j] - w[0*(NX+NU)+NX+j];
        // for ( int_t j = 0; j < NBU; j++ ) pub[0][j] = u_max[j] - w[0*(NX+NU)+NX+j];

        // dgemv_n_3l(NX, NX, pA[0], NX, x0, pb[0]);

        for (int_t j = 0; j < NX; j++) {
            pq[N][j] = Q[j*(NX+1)]*(w[N*(NX+NU)+j]-xref[j]);
        }

        // void *mem = 0;

        status = ocp_qp_hpmpc(&qp_in, &qp_out, &hpmpc_args, mem, workspace);

        // int status = 0;
        // printf("hpmpc_status=%i\n", status);

        // }
        // for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
        // shift_states(w, x_end, N);
        // shift_controls(w, u_end, N);
        timings = acados_toc(&timer);
        sum_timings+=timings;
        if (timings < min_timings) min_timings = timings;
      }
      // done with timings, store closed loop solution and
      // update primal and dual variables

      if (status == 1) printf("status = ACADOS_MAXITER\n");

      if (status == 2) printf("status = ACADOS_MINSTEP\n");

      // there is no x0 in the first stage
      for (int_t j = 0; j < NU; j++) w[0*(NX+NU)+NX+j] += qp_out.u[0][j];
      // for (int_t j = 0; j < NX; j++) pi_n[0*NX+j] = qp_out.pi[0][j];
      for (int_t j = 0; j < 2*(NBX+NBU); j++) lam_n[0*2*(NBX+NBU)+j] = qp_out.lam[0][j] + GAMMA;
      for (int_t j = 0; j < 2*(NBX+NBU); j++) t_n[0*2*(NBX+NBU)+j] = qp_out.t[0][j] + GAMMA;

      for (int_t i = 1; i < N; i++) {
          for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] += qp_out.x[i][j];
          for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] += qp_out.u[i][j];
          // for (int_t j = 0; j < NX; j++) pi_n[0*NX+j] = qp_out.pi[0][j];
          for (int_t j = 0; j < 2*(NBX+NBU); j++) lam_n[i*2*(NBX+NBU)+j] = qp_out.lam[i][j]+ GAMMA;
          for (int_t j = 0; j < 2*(NBX+NBU); j++) t_n[i*2*(NBX+NBU)+j] = qp_out.t[i][j]+ GAMMA;
      }
      for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] += qp_out.x[N][j];
      // for (int_t j = 0; j < NX; j++) pi_n[0*NX+j] = qp_out.pi[0][j];
      for (int_t j = 0; j < 2*NBX; j++) lam_n[N*2*(NBX+NBU)+j] = qp_out.lam[N][j]+ GAMMA;
      for (int_t j = 0; j < 2*NBX; j++) t_n[N*2*(NBX+NBU)+j] = qp_out.t[N][j]+ GAMMA;

      // for (int_t j = 0; j < NX; j++) w_cl[sim_iter*(NX+NU) + j] = w[j];
      // for (int_t j = 0; j < NU; j++) w_cl[sim_iter*(NX+NU) + NX + j] = w[j+NX];

      // update initial condition
      for (int_t j = 0; j < NX; j++) w[j] = w[(NX+NU) + j];
    }
    #ifdef DEBUG
     print_states_controls(&w[0], N);
    #endif  // DEBUG

    #ifdef PLOT_RESULTS
    plot_states_controls(w, T);
    #endif  // PLOT_RESULTS

    int ip_iter = hpmpc_args.out_iter;
    printf("Solved in %d iterations.\n", ip_iter);
    printf("Average of %.3f ms per RTI.\n", 1e3*sum_timings/NREP);
    printf("Minimum of %.3f ms per RTI.\n", 1e3*min_timings);
    free(workspace);
    ocp_qp_hpmpc_free_memory(mem);
    sim_erk_destroy(erk_work);
    free(sim_in.x);
    free(sim_in.u);
    free(sim_in.S_forw);
    free(sim_out.xn);
    free(sim_out.S_forw);

    // UGLY
    free(px0[0]);
    free(plb[0]);
    free(pub[0]);
    for (int_t i = 0; i < N; i++) {
        free(pA[i]);
        free(pB[i]);
        free(pb[i]);
        free(pS[i]);
        free(pq[i]);
        free(pr[i]);
        free(px[i]);
        free(pu[i]);
        free(plb[i+1]);
        free(pub[i+1]);
    }
    free(pq[N]);
    free(px[N]);

    free(A0);

    free(C0);
    free(D0);
    free(lg0);
    free(ug0);

    free(C);
    free(D);
    free(lg);
    free(ug);

    free(CN);
    free(DN);
    free(lgN);  // force all states to 0 at the last stage
    free(ugN);  // force all states to 0 at the last stage

    ii = 0;
    for (ii = 0; ii < N; ii++) {
        free(ppi[ii]);
        free(plam[ii]);
        free(pt[ii]);

        free(lam_in[ii]);
        free(t_in[ii]);
        free(ux_in[ii]);
    }

    free(plam[N]);
    free(pt[N]);

    free(lam_in[N]);
    free(t_in[N]);
    free(ux_in[N]);

    return 0;
}

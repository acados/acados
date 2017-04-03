#include "freertos/FreeRTOS.h"
#include "esp_wifi.h"
#include "esp_system.h"
#include "esp_event.h"
#include "esp_event_loop.h"
#include "nvs_flash.h"
#include "driver/gpio.h"

/* Begin acados code */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
// flush denormals to zero
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) ||  \
    defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || \
    defined(TARGET_AMD_SSE3)
#include <xmmintrin.h>  // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

#include "aux_d.h"
#include "blasfeo_d_aux.h"

#include "ocp_qp_common.h"
#include "ocp_qp_hpmpc.h"
#include "tools.h"
#include "types.h"

// define IP solver arguments && number of repetitions
#define NREP 10
#define MAXITER 10
#define TOL 1e-8
#define MINSTEP 1e-8

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
    for(ii = 0; ii < pp; ii++) T[ii * (pp + 1)] = -2;
    for(ii = 0; ii < pp - 1; ii++) T[ii * (pp + 1) + 1] = 1;
    for(ii = 1; ii < pp; ii++) T[ii * (pp + 1) - 1] = 1;

    double *Z;
    d_zeros(&Z, pp, pp);
    double *I;
    d_zeros(&I, pp, pp);
    for(ii = 0; ii < pp; ii++) I[ii * (pp + 1)] = 1.0;   // = eye(pp);
    double *Ac;
    d_zeros(&Ac, nx, nx);
    dmcopy(pp, pp, Z, pp, Ac, nx);
    dmcopy(pp, pp, T, pp, Ac + pp, nx);
    dmcopy(pp, pp, I, pp, Ac + pp * nx, nx);
    dmcopy(pp, pp, Z, pp, Ac + pp * (nx + 1), nx);
    free(T);
    free(Z);
    free(I);

    d_zeros(&I, nu, nu);
    for(ii = 0; ii < nu; ii++) I[ii * (nu + 1)] = 1.0;   // I = eye(nu);
    double *Bc;
    d_zeros(&Bc, nx, nu);
    dmcopy(nu, nu, I, nu, Bc + pp, nx);
    free(I);

    /************************************************
    * compute the discrete time system
    ************************************************/

    double *bb;
    d_zeros(&bb, nx, 1);
    dmcopy(nx, 1, bb, nx, b, nx);

    dmcopy(nx, nx, Ac, nx, A, nx);
    dscal_3l(nx2, Ts, A);
    expm(nx, A);

    d_zeros(&T, nx, nx);
    d_zeros(&I, nx, nx);
    for(ii = 0; ii < nx; ii++) I[ii * (nx + 1)] = 1.0;   // I = eye(nx);
    dmcopy(nx, nx, A, nx, T, nx);
    daxpy_3l(nx2, -1.0, I, T);
    dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);

    int *ipiv = (int *)malloc(nx * sizeof(int));
    dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
    free(ipiv);

    free(Ac);
    free(Bc);
    free(bb);

    /************************************************
    * initial state
    ************************************************/

    if(nx == 4) {
        x0[0] = 5;
        x0[1] = 10;
        x0[2] = 15;
        x0[3] = 20;
    } else {
        int jj;
        for(jj = 0; jj < nx; jj++) x0[jj] = 1;
    }
}

/* End acados code */

esp_err_t event_handler(void *ctx, system_event_t *event)
{
    return ESP_OK;
}

void main_memory_task(void *pv)
{
  // result on board Nano32 running ESP32:
  //  More stack size == less heap size for malloc(). Their sum is about 180 kB.

  unsigned long occupy_heap_size=1024;

  // void *dump_heap_size = malloc(occupy_heap_size); // for debug 115111: already big
  while(1) {
    int *dump_heap_size = malloc(occupy_heap_size);
    if(dump_heap_size == NULL) {
      printf("Pointer %p, failed to allocate %lu bytes.\n", dump_heap_size, occupy_heap_size);
      break;
    }
    else {
    printf("Pointer %p, Allocated %lu bytes.\n", dump_heap_size, occupy_heap_size);
    free(dump_heap_size);
    occupy_heap_size += 1024;
    }
  }
}

void main_task(void *pv)
{
  int occupy_heap_size = 159000; // for this problem, HPMPC with partial condensing need about 158xxx bytes for work space
  void *dump_heap_size = malloc(occupy_heap_size); // for debug 115111: already big
  printf("\n Allocated %d bytes, at pointer %p",occupy_heap_size,dump_heap_size);
    int loopnumber = 0; // for debug

    /* Begin acados code */
    printf("\n");
    printf("\n");
    printf("\n");
    printf(
        " HPMPC -- Library for High-Performance implementation of solvers for "
        "MPC.\n");
    printf(
        " Copyright (C) 2014-2015 by Technical University of Denmark. All "
        "rights reserved.\n");
    printf("\n");
    printf(" HPMPC is distributed in the hope that it will be useful,\n");
    printf(" but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf(" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
    printf(" See the GNU Lesser General Public License for more details.\n");
    printf("\n");
    printf("\n");
    printf("\n");

#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) ||  \
        defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || \
        defined(TARGET_AMD_SSE3)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);  // flush to zero subnormals !!!
    // works only with one thread
    // !!!
#endif

    int ii, jj;

    int rep, nrep = NREP;

    int nx = 8;  // number of states (it has to be even for the mass-spring
    // system test problem)
    int nu = 3;  // number of inputs (controllers) (it has to be at least 1 and
    // at most nx/2 for the mass-spring system test problem)
    int N = 15;   // horizon length
    int nb = 11;  // number of box constrained inputs and states
    int ng = 0;  // 4;  // number of general constraints
    int ngN = 4;  // 4;  // number of general constraints at the last stage

    int N2 = 3;   // horizon length

    int nbu = nu < nb ? nu : nb;
    int nbx = nb - nu > 0 ? nb - nu : 0;

    // stage-wise variant size
    int nxx[N + 1];
    nxx[0] = 0;
    for(ii = 1; ii <= N; ii++) nxx[ii] = nx;

    int nuu[N + 1];
    for(ii = 0; ii < N; ii++) nuu[ii] = nu;
    nuu[N] = 0;

    int nbb[N + 1];
    nbb[0] = nbu;
    for(ii = 1; ii < N; ii++) nbb[ii] = nb;
    nbb[N] = nbx;

    int ngg[N + 1];
    for(ii = 0; ii < N; ii++) ngg[ii] = ng;
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
    printf(
        " IP method parameters: predictor-corrector IP, double precision, %d "
        "maximum iterations, %5.1e exit tolerance in duality measure.\n",
        MAXITER, TOL);
    printf("\n");
#if defined(TARGET_X64_AVX2)
    printf(" HPMPC built for the AVX2 architecture\n");
#endif
#if defined(TARGET_X64_AVX)
    printf(" HPMPC built for the AVX architecture\n");
#endif
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

    for(jj = 0; jj < nx; jj++) b[jj] = 0.1;

    for(jj = 0; jj < nx; jj++) x0[jj] = 0;
    x0[0] = 2.5;
    x0[1] = 2.5;

    //    d_print_mat(nx, nx, A, nx);
    //    d_print_mat(nx, nu, B, nx);
    //    d_print_mat(nx, 1, b, nx);
    //    d_print_mat(nx, 1, x0, nx);

    // compute b0 = b + A*x0
    double *b0;
    d_zeros(&b0, nx, 1);
    dcopy_3l(nx, b, 1, b0, 1);
    dgemv_n_3l(nx, nx, A, nx, x0, b0);
    //    d_print_mat(nx, 1, b, nx);
    //    d_print_mat(nx, 1, b0, nx);

    // then A0 is a matrix of size 0x0
    double *A0;
    d_zeros(&A0, 0, 0);

    /************************************************
    * box constraints
    ************************************************/

    int *idxb0;
    int_zeros(&idxb0, nbb[0], 1);
    double *lb0;
    d_zeros(&lb0, nbb[0], 1);
    double *ub0;
    d_zeros(&ub0, nbb[0], 1);
    for(jj = 0; jj < nbu; jj++) {
        lb0[jj] = -0.5;  //   umin
        ub0[jj] = 0.5;   //   umax
        idxb0[jj] = nxx[0] + jj;
    }
    //    int_print_mat(nbb[0], 1, idxb0, nbb[0]);

    int *idxb1;
    int_zeros(&idxb1, nbb[1], 1);
    double *lb1;
    d_zeros(&lb1, nbb[1], 1);
    double *ub1;
    d_zeros(&ub1, nbb[1], 1);
    for(jj = 0; jj < nbx; jj++) {
        lb1[jj] = -4.0;  //   xmin
        ub1[jj] = 4.0;   //   xmax
        idxb1[jj] = jj;
    }
    for(; jj < nb; jj++) {
        lb1[jj] = -0.5;  //   umin
        ub1[jj] = 0.5;   //   umax
        idxb1[jj] = jj;
    }
    //    int_print_mat(nbb[1], 1, idxb1, nbb[1]);

    int *idxbN;
    int_zeros(&idxbN, nbb[N], 1);
    double *lbN;
    d_zeros(&lbN, nbb[N], 1);
    double *ubN;
    d_zeros(&ubN, nbb[N], 1);
    for(jj = 0; jj < nbx; jj++) {
        lbN[jj] = -4.0;  //   umin
        ubN[jj] = 4.0;   //   umax
        idxbN[jj] = jj;
    }
    //    int_print_mat(nbb[N], 1, idxb1, nbb[N]);

    /************************************************
    * general constraints
    ************************************************/

    double *C;
    d_zeros(&C, ng, nx);
    double *D;
    d_zeros(&D, ng, nu);
    double *lg;
    d_zeros(&lg, ng, 1);
    double *ug;
    d_zeros(&ug, ng, 1);

    double *CN;
    d_zeros(&CN, ngN, nx);
    for(ii = 0; ii < ngN; ii++) CN[ii * (ngN + 1)] = 1.0;
    //    d_print_mat(ngN, nx, CN, ngN);
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

    /************************************************
    * cost function
    ************************************************/

    double *Q;
    d_zeros(&Q, nx, nx);
    for(ii = 0; ii < nx; ii++) Q[ii * (nx + 1)] = 1.0;

    double *R;
    d_zeros(&R, nu, nu);
    for(ii = 0; ii < nu; ii++) R[ii * (nu + 1)] = 2.0;

    double *S;
    d_zeros(&S, nu, nx);

    double *q;
    d_zeros(&q, nx, 1);
    for(ii = 0; ii < nx; ii++) q[ii] = 0.1;

    double *r;
    d_zeros(&r, nu, 1);
    for(ii = 0; ii < nu; ii++) r[ii] = 0.2;

    // Q0 and q0 are matrices of size 0
    double *Q0;
    d_zeros(&Q0, 0, 0);
    double *q0;
    d_zeros(&q0, 0, 1);

    // compute r0 = r + S*x0
    double *r0;
    d_zeros(&r0, nu, 1);
    dcopy_3l(nu, r, 1, r0, 1);
    dgemv_n_3l(nu, nx, S, nu, x0, r0);

    // then S0 is a matrix of size nux0
    double *S0;
    d_zeros(&S0, nu, 0);

    /* End acados code */

    gpio_set_direction(GPIO_NUM_16, GPIO_MODE_OUTPUT);
    int level = 0;
    while(true) {
        gpio_set_level(GPIO_NUM_16, level);
        level = !level;
        vTaskDelay(100 / portTICK_PERIOD_MS);
        printf("\n\n New loop %d\n\n", ++loopnumber);

        /* Begin acados code */

        /************************************************
        * problems data
        ************************************************/

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

        hA[0] = A0;
        hB[0] = B;
        hb[0] = b0;
        hQ[0] = Q0;
        hS[0] = S0;
        hR[0] = R;
        hq[0] = q0;
        hr[0] = r0;
        hlb[0] = lb0;
        hub[0] = ub0;
        hidxb[0] = idxb0;
        hC[0] = C;
        hD[0] = D;
        hlg[0] = lg;
        hug[0] = ug;
        for(ii = 1; ii < N; ii++) {
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
        double *hu[N];
        double *hpi[N];
        double *hlam[N + 1];

        for(ii = 0; ii < N; ii++) {
            d_zeros(&hx[ii], nxx[ii], 1);
            d_zeros(&hu[ii], nuu[ii], 1);
            d_zeros(&hpi[ii], nxx[ii + 1], 1);
            d_zeros(&hlam[ii], 2 * nbb[ii] + 2 * nbb[ii], 1); // Andrea: why do we have 4*nb here?
        }
        d_zeros(&hx[N], nxx[N], 1);
        d_zeros(&hlam[N], 2 * nbb[N] + 2 * nbb[N], 1);

        /************************************************
        * create the in and out struct
        ************************************************/

        ocp_qp_in qp_in;
        qp_in.N = N;
        qp_in.nx = (const int *) nxx;
        qp_in.nu = (const int *) nuu;
        qp_in.nb = (const int *) nbb;
        qp_in.nc = (const int *) ngg;
        qp_in.A = (const double **) hA;
        qp_in.B = (const double **) hB;
        qp_in.b = (const double **) hb;
        qp_in.Q = (const double **) hQ;
        qp_in.S = (const double **) hS;
        qp_in.R = (const double **) hR;
        qp_in.q = (const double **) hq;
        qp_in.r = (const double **) hr;
        qp_in.idxb = (const int **) hidxb;
        qp_in.lb = (const double **) hlb;
        qp_in.ub = (const double **) hub;
        qp_in.Cx = (const double **) hC;
        qp_in.Cu = (const double **) hD;
        qp_in.lc = (const double **) hlg;
        qp_in.uc = (const double **) hug;

        ocp_qp_out qp_out;
        qp_out.x = hx;
        qp_out.u = hu;
        qp_out.pi = hpi;
        qp_out.lam = hlam;

        /************************************************
        * solver arguments (fully sparse)
        ************************************************/

        // solver arguments
        ocp_qp_hpmpc_args hpmpc_args;
        hpmpc_args.tol = TOL;
        hpmpc_args.max_iter = MAXITER;
        //  hpmpc_args.min_step = MINSTEP;
        hpmpc_args.mu0 = 0.0;
        //  hpmpc_args.sigma_min = 1e-3;
        hpmpc_args.warm_start = 0;
        hpmpc_args.N2 = N;
        double inf_norm_res[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        hpmpc_args.inf_norm_res = &inf_norm_res[0];

        /************************************************
        * work space (fully sparse)
        ************************************************/
        int return_value;
        struct timeval tv0, tv1;
if (0) {
        int work_space_size =
            ocp_qp_hpmpc_workspace_size_bytes(N, nxx, nuu, nbb, ngg, hidxb, &hpmpc_args);
        printf("\nwork space size: %d bytes\n", work_space_size);
        printf("Free heap size: %d\n",esp_get_free_heap_size()); // for debug
        free(dump_heap_size); // free the contiguous block in the heap created before
        void *workspace = malloc(work_space_size);
        printf("\nPointer to workspace: workspace=%p\n", workspace);

        /************************************************
        * call the solver (fully sparse)
        ************************************************/

        // int return_value;

        // struct timeval tv0, tv1;
        gettimeofday(&tv0, NULL);  // stop
        printf("gettimeofday, line 544\n"); // debug

        for(rep = 0; rep < nrep; rep++) {
            printf("Iteration %d \n", rep); // debug

            // call the QP OCP solver
            return_value = ocp_qp_hpmpc(&qp_in, &qp_out, &hpmpc_args, workspace);
        }

        gettimeofday(&tv1, NULL);  // stop
        printf("gettimeofday, line 554\n"); // debug

        if(return_value == ACADOS_SUCCESS)
            printf("\nACADOS status: solution found in %d iterations\n", hpmpc_args.out_iter);

        if(return_value == ACADOS_MAXITER)
            printf("\nACADOS status: maximum number of iterations reached\n");

        if(return_value == ACADOS_MINSTEP)
            printf("\nACADOS status: below minimum step size length\n");

        // Debug
        printf("\nu = \n");
        for (ii = 0; ii < N; ii++) d_print_mat(1, nuu[ii], hu[ii], 1);

        printf("\nx = \n");
        for (ii = 0; ii <= N; ii++) d_print_mat(1, nxx[ii], hx[ii], 1);

        double time = (tv1.tv_sec - tv0.tv_sec) / (nrep + 0.0) +
                      (tv1.tv_usec - tv0.tv_usec) / (nrep * 1e6);

        printf("\n");
        printf(" inf norm res: %e, %e, %e, %e, %e\n", inf_norm_res[0], inf_norm_res[1], \
               inf_norm_res[2], inf_norm_res[3], inf_norm_res[4]);
        printf("\n");
        printf(" Average solution time over %d runs: %5.2e seconds\n", nrep, time);
        printf("\n\n");
        // free allocated memory for work space
        free(workspace);
}
        /************************************************
        * solver arguments (partial condensing)
        ************************************************/
        // void *dump_heap_size = malloc(occupy_heap_size); // for debug 115111: already big
        // printf("\n Allocated %d bytes, at pointer %p",occupy_heap_size,dump_heap_size);

        // solver arguments
        hpmpc_args.N2 = N2;


        /************************************************
        * work space (partial condensing)
        ************************************************/

        int work_space_size_part_cond =
            ocp_qp_hpmpc_workspace_size_bytes(N, nxx, nuu, nbb, ngg, hidxb, &hpmpc_args);
        printf("\nwork space size: %d bytes\n", work_space_size_part_cond);
        printf("Free heap size: %d\n",esp_get_free_heap_size()); // for debug 231296
        free(dump_heap_size); // free the contiguous block in the heap created before
        void *workspace_part_cond = malloc(work_space_size_part_cond);
        // void *workspace_part_cond = malloc(123456); // for debug 115111: already big
        // printf("Work space to alloc: 123456\n"); // for debug
        printf("\nPointer to workspace: workspace_part_cond=%p\n", workspace_part_cond);
        printf("Free heap size left: %d\n",esp_get_free_heap_size()); // for debug

        /************************************************
        * call the solver (partial condensing)
        ************************************************/

        gettimeofday(&tv0, NULL);  // stop

        for(rep = 0; rep < nrep; rep++) {
            // call the QP OCP solver
            return_value = ocp_qp_hpmpc(&qp_in, &qp_out, &hpmpc_args, workspace_part_cond);
        }

        gettimeofday(&tv1, NULL);  // stop

        if(return_value == ACADOS_SUCCESS)
            printf("\nACADOS status: solution found in %d iterations\n", hpmpc_args.out_iter);

        if(return_value == ACADOS_MAXITER)
            printf("\nACADOS status: maximum number of iterations reached\n");

        if(return_value == ACADOS_MINSTEP)
            printf("\nACADOS status: below minimum step size length\n");

        printf("\nu = \n");
        for(ii = 0; ii < N; ii++) d_print_mat(1, nuu[ii], hu[ii], 1);

        printf("\nx = \n");
        for(ii = 0; ii <= N; ii++) d_print_mat(1, nxx[ii], hx[ii], 1);

        double time_part_cond = (tv1.tv_sec - tv0.tv_sec) / (nrep + 0.0) +
                                (tv1.tv_usec - tv0.tv_usec) / (nrep * 1e6);

        printf("\n");
        printf(" inf norm res: %e, %e, %e, %e, %e\n", inf_norm_res[0], inf_norm_res[1], \
               inf_norm_res[2], inf_norm_res[3], inf_norm_res[4]);
        printf("\n");
        printf(" Average solution time over %d runs (part cond): %5.2e seconds\n", nrep, \
               time_part_cond);
        printf("\n\n");

        /************************************************
        * free memory
        ************************************************/

        // d_free(A);
        // d_free(B);
        // d_free(b);
        // d_free(x0);
        // d_free(A0);
        // d_free(b0);
        // d_free(Q);
        // d_free(S);
        // d_free(R);
        // d_free(q);
        // d_free(r);
        // d_free(Q0);
        // d_free(S0);
        // d_free(q0);
        // d_free(r0);
        // int_free(idxb0);
        // d_free(lb0);
        // d_free(ub0);
        // int_free(idxb1);
        // d_free(lb1);
        // d_free(ub1);
        // int_free(idxbN);
        // d_free(lbN);
        // d_free(ubN);
        // d_free(C);
        // d_free(D);
        // d_free(lg);
        // d_free(ug);
        // d_free(CN);
        // d_free(lgN);
        // d_free(ugN);

        // for(ii = 0; ii < N; ii++) {
        //     d_free(hx[ii]);
        //     d_free(hu[ii]);
        // }
        // d_free(hx[N]);

        // free(workspace);
        // free(workspace_part_cond);

        /* End acados code */

    }
}
void app_main(void)
{
    nvs_flash_init();
    // tcpip_adapter_init();
    // ESP_ERROR_CHECK( esp_event_loop_init(event_handler, NULL) );
    // wifi_init_config_t cfg = WIFI_INIT_CONFIG_DEFAULT();
    // ESP_ERROR_CHECK( esp_wifi_init(&cfg) );
    // ESP_ERROR_CHECK( esp_wifi_set_storage(WIFI_STORAGE_RAM) );
    // ESP_ERROR_CHECK( esp_wifi_set_mode(WIFI_MODE_STA) );
    // wifi_config_t sta_config = {
    //     .sta = {
    //         .ssid = "access_point_name",
    //         .password = "password",
    //         .bssid_set = false
    //     }
    // };
    // ESP_ERROR_CHECK( esp_wifi_set_config(WIFI_IF_STA, &sta_config) );
    // ESP_ERROR_CHECK( esp_wifi_start() );
    // ESP_ERROR_CHECK( esp_wifi_connect() );
    xTaskCreate(main_task, "main", 10*1024, NULL, 5, NULL);
    // xTaskCreate(main_memory_task, "main_memory", 10*1024, NULL, 5, NULL);
}

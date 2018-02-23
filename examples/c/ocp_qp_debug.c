
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
// #include <unistd.h>  // NOTE(dimitris): to read current directory

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/print.h"
#include "acados_c/ocp_qp.h"

// needed for casting args
#include "acados/ocp_qp/ocp_qp_sparse_solver.h"
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_hpipm.h"
#include "acados/dense_qp/dense_qp_qpoases.h"

#ifdef ACADOS_WITH_HPMPC
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#endif

#ifdef ACADOS_WITH_QPDUNES
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#endif

#ifdef ACADOS_WITH_QORE
#include "acados/dense_qp/dense_qp_qore.h"
#endif

// TODOS!
// STORE OPTIMAL SOLUTION FROM MATLAB (and DEFINE with SOLVER)
// RUN CLOSED LOOP SIMULATION AND PRINT NITER, CPU_TIME, ERROR
// NREP TO TAKE MINIMUM TIME

void load_ptr(void *lib, char *data_string, void **ptr)
{
    *ptr = dlsym(lib, data_string);
    if (ptr == NULL)
    {
        printf("dlsym failed: %s\n", dlerror());
        exit(1);
    }
}

int main() {

    int n_problems = 5;

    /************************************************
    * load dynamic library
    ************************************************/

    char str[256];

    // TODO(dimitris): currently assuming we run it from build dir
    snprintf(str, sizeof(str), "../examples/c/ocp_qp_bugs/ocp_qp_data.so");

    void *lib = dlopen(str, RTLD_NOW);
    if (lib == NULL) {
        printf("dlopen failed: %s\n", dlerror());
        exit(1);
    }

    /************************************************
    * set up dims struct
    ************************************************/

    ocp_qp_dims dims;

    int N, *N_ptr;

    snprintf(str, sizeof(str), "N");
    load_ptr(lib, str, (void **)&N_ptr);

    N = *N_ptr;
    dims.N = N;

    snprintf(str, sizeof(str), "nx");
    load_ptr(lib, str, (void **)&dims.nx);

    snprintf(str, sizeof(str), "nu");
    load_ptr(lib, str, (void **)&dims.nu);

    snprintf(str, sizeof(str), "nb");
    load_ptr(lib, str, (void **)&dims.nb);

    snprintf(str, sizeof(str), "nbu");
    load_ptr(lib, str, (void **)&dims.nbu);

    snprintf(str, sizeof(str), "nbx");
    load_ptr(lib, str, (void **)&dims.nbx);

    snprintf(str, sizeof(str), "ng");
    load_ptr(lib, str, (void **)&dims.ng);

    snprintf(str, sizeof(str), "ns");
    load_ptr(lib, str, (void **)&dims.ns);

    int **hidxb = malloc((N+1)*sizeof(int *));

    int *idxb;

    snprintf(str, sizeof(str), "idxb");
    load_ptr(lib, str, (void **)&idxb);

    int sum_idxb = 0;

    for (int ii = 0; ii < N+1; ii++)
    {
        hidxb[ii] = &idxb[sum_idxb];
        sum_idxb += dims.nb[ii];
    }

    /************************************************
    * set up solver and input/output
    ************************************************/

    ocp_qp_solver_plan plan;
    plan.qp_solver = FULL_CONDENSING_QORE;

    void *args = ocp_qp_create_args(&plan, &dims);

    int N2 = dims.N;

    ocp_qp_full_condensing_solver_args *fcond_solver_args = NULL;
    ocp_qp_sparse_solver_args *pcond_solver_args = NULL;
    pcond_solver_args++;
    fcond_solver_args++;

    switch (plan.qp_solver)
    {
    case PARTIAL_CONDENSING_HPIPM:
        printf("\nPartial condensing + HPIPM (N2 = %d):\n\n", N2);
        pcond_solver_args = (ocp_qp_sparse_solver_args *)args;
        ocp_qp_partial_condensing_args *hpipm_pcond_args = (ocp_qp_partial_condensing_args *)pcond_solver_args->pcond_args;
        ocp_qp_hpipm_args *hpipm_solver_args = (ocp_qp_hpipm_args *)pcond_solver_args->solver_args;

        hpipm_pcond_args->N2 = N2;
        hpipm_solver_args->hpipm_args->iter_max = 1000;
        hpipm_solver_args->hpipm_args->warm_start = 1;  // // TODO(dimitris): ONLY WORKS ONLY WITH WARM_START!
        // hpipm_solver_args->hpipm_args->mu0 = 1e6;
        // hpipm_solver_args->hpipm_args->tol = 1e-12;

        break;
    case PARTIAL_CONDENSING_HPMPC:
#ifdef ACADOS_WITH_HPMPC
        printf("\nPartial condensing + HPMPC (N2 = %d):\n\n", N2);
        pcond_solver_args = (ocp_qp_sparse_solver_args *)args;
        ocp_qp_partial_condensing_args *hpmpc_pcond_args = (ocp_qp_partial_condensing_args *)pcond_solver_args->pcond_args;
        ocp_qp_hpmpc_args *hpmpc_solver_args = (ocp_qp_hpmpc_args *)pcond_solver_args->solver_args;

        hpmpc_pcond_args->N2 = N2;

        hpmpc_solver_args->max_iter = 1000;
        hpmpc_solver_args->tol = 1e-12;
#endif
        break;
    case PARTIAL_CONDENSING_QPDUNES:
#ifdef ACADOS_WITH_QPDUNES
        printf("\nPartial condensing + qpDUNES (N2 = %d):\n\n", N2);
        pcond_solver_args = (ocp_qp_sparse_solver_args *)args;
        ocp_qp_partial_condensing_args *qpdunes_pcond_args = (ocp_qp_partial_condensing_args *)pcond_solver_args->pcond_args;
        ocp_qp_qpdunes_args *qpdunes_solver_args = (ocp_qp_qpdunes_args *)pcond_solver_args->solver_args;

        qpdunes_pcond_args->N2 = N2;  // NOTE(dimitris): only change N2 above, not here!

        qpdunes_solver_args->warmstart = 1;

        if (N2 == dims.N)
        {
            qpdunes_solver_args->stageQpSolver = QPDUNES_WITH_CLIPPING;
            qpdunes_solver_args->options.lsType = QPDUNES_LS_ACCELERATED_GRADIENT_BISECTION_LS;
            // NOTE(dimitris): these two options should always change together
        } else
        {
            qpdunes_solver_args->stageQpSolver = QPDUNES_WITH_QPOASES;
            qpdunes_solver_args->options.lsType = QPDUNES_LS_HOMOTOPY_GRID_SEARCH;
        }
#endif
        break;
    case FULL_CONDENSING_HPIPM:
        printf("\nFull condensing + HPIPM:\n\n");
        fcond_solver_args = (ocp_qp_full_condensing_solver_args *)args;
        dense_qp_hpipm_args *dense_hpipm_solver_args = (dense_qp_hpipm_args *)fcond_solver_args->solver_args;

        dense_hpipm_solver_args->hpipm_args->iter_max = 1000;
        dense_hpipm_solver_args->hpipm_args->warm_start = 0;

        break;
    case FULL_CONDENSING_QORE:
#ifdef ACADOS_WITH_QORE
        printf("\nFull condensing + QORE:\n\n");
        fcond_solver_args = (ocp_qp_full_condensing_solver_args *)args;
        dense_qp_qore_args *qore_solver_args = (dense_qp_qore_args *)fcond_solver_args->solver_args;

        qore_solver_args->max_iter = 1000;
        qore_solver_args->warm_start = 0;  // TODO(dimitris): ONLY WORKS COLD STARTED

        if (qore_solver_args->warm_start)
            qore_solver_args->warm_strategy = 0;

        break;
#endif
    case FULL_CONDENSING_QPOASES:
        printf("\nFull condensing + QPOASES:\n\n");
        fcond_solver_args = (ocp_qp_full_condensing_solver_args *)args;
        dense_qp_qpoases_args *qpoases_solver_args = (dense_qp_qpoases_args *)fcond_solver_args->solver_args;

        qpoases_solver_args->warm_start = 1;

        break;
    case PARTIAL_CONDENSING_OOQP:
        break;
    }

    ocp_qp_in *qp_in = create_ocp_qp_in(&dims);

    ocp_qp_out *qp_out = create_ocp_qp_out(&dims);

    ocp_qp_solver *qp_solver = ocp_qp_create(&plan, &dims, args);

    /************************************************
    * define pointers to be used in closed loop
    ************************************************/

    // dynamics
    double **hA = malloc(N*sizeof(double *));
    double **hB = malloc(N*sizeof(double *));
    double **hb = malloc(N*sizeof(double *));

    double *Av, *Bv, *b;

    int sum_A, sum_B, sum_b;

    // objective
    double **hQ = malloc((N+1)*sizeof(double *));
    double **hR = malloc((N+1)*sizeof(double *));
    double **hS = malloc((N+1)*sizeof(double *));
    double **hq = malloc((N+1)*sizeof(double *));
    double **hr = malloc((N+1)*sizeof(double *));

    double *Qv, *Rv, *Sv, *q, *r;

    int sum_Q, sum_R, sum_S, sum_q, sum_r;

    // constraints
    double **hlb = malloc((N+1)*sizeof(double *));
    double **hub = malloc((N+1)*sizeof(double *));
    double **hC = malloc((N+1)*sizeof(double *));
    double **hD = malloc((N+1)*sizeof(double *));
    double **hlg = malloc((N+1)*sizeof(double *));
    double **hug = malloc((N+1)*sizeof(double *));

    double *lb, *ub;

    int sum_nb = 0;

    // info
    ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;

    for (int indx = 0; indx < n_problems; indx++)
    {

        /************************************************
        * set up dynamics
        ************************************************/

        snprintf(str, sizeof(str), "Av_%d", indx);
        load_ptr(lib, str, (void **)&Av);

        snprintf(str, sizeof(str), "Bv_%d", indx);
        load_ptr(lib, str, (void **)&Bv);

        snprintf(str, sizeof(str), "b_%d", indx);
        load_ptr(lib, str, (void **)&b);

        sum_A = 0;
        sum_B = 0;
        sum_b = 0;

        for (int ii = 0; ii < N; ii++)
        {
            // TODO(dimitris): test for varying dimensions
            hA[ii] = &Av[sum_A];
            sum_A += dims.nx[ii+1]*dims.nx[ii];

            hB[ii] = &Bv[sum_B];
            sum_B += dims.nx[ii+1]*dims.nu[ii];

            hb[ii] = &b[sum_b];
            sum_b += dims.nx[ii+1];
        }

        /************************************************
        * set up objective
        ************************************************/

        snprintf(str, sizeof(str), "Qv_%d", indx);
        load_ptr(lib, str, (void **)&Qv);

        snprintf(str, sizeof(str), "Rv_%d", indx);
        load_ptr(lib, str, (void **)&Rv);

        snprintf(str, sizeof(str), "Sv_%d", indx);
        load_ptr(lib, str, (void **)&Sv);

        snprintf(str, sizeof(str), "q_%d", indx);
        load_ptr(lib, str, (void **)&q);

        snprintf(str, sizeof(str), "r_%d", indx);
        load_ptr(lib, str, (void **)&r);

        sum_Q = 0;
        sum_R = 0;
        sum_S = 0;
        sum_q = 0;
        sum_r = 0;

        for (int ii = 0; ii < N+1; ii++)
        {
            hQ[ii] = &Qv[sum_Q];
            sum_Q += dims.nx[ii]*dims.nx[ii];

            hR[ii] = &Rv[sum_R];
            sum_R += dims.nu[ii]*dims.nu[ii];

            hS[ii] = &Sv[sum_S];
            sum_S += dims.nx[ii]*dims.nu[ii];

            hq[ii] = &q[sum_q];
            sum_q += dims.nx[ii];

            hr[ii] = &r[sum_r];
            sum_r += dims.nu[ii];
        }

        /************************************************
        * set up constraints
        ************************************************/

        snprintf(str, sizeof(str), "lb_%d", indx);
        load_ptr(lib, str, (void **)&lb);

        snprintf(str, sizeof(str), "ub_%d", indx);
        load_ptr(lib, str, (void **)&ub);

        sum_nb = 0;

        for (int ii = 0; ii < N+1; ii++)
        {
            hlb[ii] = &lb[sum_nb];
            hub[ii] = &ub[sum_nb];
            sum_nb += dims.nb[ii];
        }

        /************************************************
        * convert data and solve qp
        ************************************************/

        d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, qp_in);

        // print_ocp_qp_in(qp_in);

        int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

        // print_ocp_qp_out(qp_out);

        if (acados_return != ACADOS_SUCCESS)
        {
            printf("QP SOLVER FAILED WITH FLAG %d\n", acados_return);
            return -1;
        }

        printf("\n--> problem %d solved in %d iterations\n\n", indx, info->num_iter);

    }

    /************************************************
    * free memory
    ************************************************/

    free(hA);
    free(hB);
    free(hb);

    free(hQ);
    free(hR);
    free(hS);
    free(hq);
    free(hr);

    free(hidxb);
    free(hlb);
    free(hub);
    free(hC);
    free(hD);
    free(hlg);
    free(hug);

    free(qp_in);
    free(qp_out);
    free(qp_solver);

    return 0;
}

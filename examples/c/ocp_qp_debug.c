
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <string.h>
#include <math.h>
#include <xmmintrin.h>

// #include <unistd.h>  // NOTE(dimitris): to read current directory

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/print.h"
#include "acados_c/ocp_qp.h"

#include "blasfeo_target.h"
#include "blasfeo_common.h"
#include "blasfeo_d_aux.h"

// needed for casting args
// TODO(dimitris): use options_i instead
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

#include <assert.h>

void load_ptr(void *lib, char *data_string, void **ptr)
{
    *ptr = dlsym(lib, data_string);
    if (ptr == NULL)
    {
        printf("dlsym failed: %s\n", dlerror());
        exit(1);
    }
}



void convert_strvecs_to_single_vec(int n, struct blasfeo_dvec sv[], double *v)
{
    int ind = 0;
    for (int i = 0; i < n; i++)
    {
        blasfeo_unpack_dvec(sv[i].m, &sv[i], 0, &v[ind]);
        ind += sv[i].m;
    }
}



double compare_with_acado_solution(int N, int nvars, ocp_qp_out *qp_out, double *acado_sol)
{
    double *acados_sol = malloc(nvars*sizeof(double));

    convert_strvecs_to_single_vec(N+1, qp_out->ux, acados_sol);

    double error = 0;
    double diff;

    for (int ii = 0; ii < nvars; ii++)
    {
        diff = acado_sol[ii] - acados_sol[ii];
        if (diff < 0) diff = - diff;

        if (diff > error) error = diff;
        // printf(" %2.5e\t %2.5e\n", acado_sol[ii], acados_sol[ii]);
        if isnan(acados_sol[ii])
        {
            printf("nans detected in acados solution.\n");
            exit(-1);
    }
    }

    free(acados_sol);

    return error;
}


void choose_solver(int N, char *lib_str, int *N2, int *warm_start, ocp_qp_solver_t *qp_solver);


int main() {
    // Uncomment to detect NaNs
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    int n_rep = 5;  // TODO number of runs (taking minimum time)
    int n_problems = 25;  // number of MPC problems stored in shared library

    bool auto_choose_acados_solver = true;  // choose acados solver based on lib name

    bool eliminate_x0 = false;

    char suffix[256];

    if (eliminate_x0)
        snprintf(suffix, sizeof(suffix), "_new");
    else
        snprintf(suffix, sizeof(suffix), "");

    /************************************************
    * load dynamic library
    ************************************************/

    char lib_str[256];

    char solver_in[256] = "qpDUNES_B0";
    int nmasses_in = 7;
    int warmstart_in = 1;

    int N_ins[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    int N_in;

    for (int jj = 0; jj < 10; jj++)
    {
    N_in = N_ins[jj];

    // TODO(dimitris): currently assuming we run it from build dir
    snprintf(lib_str, sizeof(lib_str),
        "../examples/c/ocp_qp_bugs/ocp_qp_data_nmasses_%d_nsteps_%d_solver_%s_warmstart_%d.so",
        nmasses_in, N_in, solver_in, warmstart_in);

    void *lib = dlopen(lib_str, RTLD_NOW);
    if (lib == NULL) {
        printf("dlopen failed: %s\n", dlerror());
        exit(1);
    }

    /************************************************
    * set up dims struct
    ************************************************/

    char str[256];

    ocp_qp_dims dims;

    int N, *N_ptr;

    snprintf(str, sizeof(str), "N");
    load_ptr(lib, str, (void **)&N_ptr);

    N = *N_ptr;
    dims.N = N;

    snprintf(str, sizeof(str), "nx%s", suffix);
    load_ptr(lib, str, (void **)&dims.nx);

    snprintf(str, sizeof(str), "nu");
    load_ptr(lib, str, (void **)&dims.nu);

    snprintf(str, sizeof(str), "nb%s", suffix);
    load_ptr(lib, str, (void **)&dims.nb);

    snprintf(str, sizeof(str), "nbu");
    load_ptr(lib, str, (void **)&dims.nbu);

    snprintf(str, sizeof(str), "nbx%s", suffix);
    load_ptr(lib, str, (void **)&dims.nbx);

    snprintf(str, sizeof(str), "ng");
    load_ptr(lib, str, (void **)&dims.ng);

    snprintf(str, sizeof(str), "ns");
    load_ptr(lib, str, (void **)&dims.ns);

    int **hidxb = malloc((N+1)*sizeof(int *));

    int *idxb;

    snprintf(str, sizeof(str), "idxb%s", suffix);
    load_ptr(lib, str, (void **)&idxb);

    int sum_idxb = 0;

    for (int ii = 0; ii < N+1; ii++)
    {
        hidxb[ii] = &idxb[sum_idxb];
        sum_idxb += dims.nb[ii];
    }

    /************************************************
    * set up solver args and input/output
    ************************************************/

    ocp_qp_solver_plan plan;
    int N2 = -1;
    int warmstart = -1;

    if (auto_choose_acados_solver == false)
    {
        // choose custom values
        N2 = dims.N;
        plan.qp_solver = FULL_CONDENSING_QPOASES;
        warmstart = 1;
    } else
    {
        // infer values from libname
        choose_solver(N, lib_str, &N2, &warmstart, &plan.qp_solver);
        if (N2 == 0) N2 = dims.N;
    }

    void *args = ocp_qp_create_args(&plan, &dims);

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
        hpipm_solver_args->hpipm_args->warm_start = warmstart;  // TODO(dimitris): ONLY WORKS ONLY WITH WARM_START!
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
        hpmpc_solver_args->warm_start = warmstart;
#endif
        break;
    case PARTIAL_CONDENSING_QPDUNES:
#ifdef ACADOS_WITH_QPDUNES
        printf("\nPartial condensing + qpDUNES (N2 = %d):\n\n", N2);
        pcond_solver_args = (ocp_qp_sparse_solver_args *)args;
        ocp_qp_partial_condensing_args *qpdunes_pcond_args = (ocp_qp_partial_condensing_args *)pcond_solver_args->pcond_args;
        ocp_qp_qpdunes_args *qpdunes_solver_args = (ocp_qp_qpdunes_args *)pcond_solver_args->solver_args;

        qpdunes_pcond_args->N2 = N2;  // NOTE(dimitris): only change N2 above, not here!

        qpdunes_solver_args->warmstart = warmstart;

        qpdunes_solver_args->options.maxIter = 1000;

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
        if (eliminate_x0)
        {
            printf("qpDUNES does not support elimination of x0, turn off flag.\n\n");
            exit(-1);
        }
#endif
        break;
    case FULL_CONDENSING_HPIPM:
        printf("\nFull condensing + HPIPM:\n\n");
        fcond_solver_args = (ocp_qp_full_condensing_solver_args *)args;
        dense_qp_hpipm_args *dense_hpipm_solver_args = (dense_qp_hpipm_args *)fcond_solver_args->solver_args;

        dense_hpipm_solver_args->hpipm_args->iter_max = 1000;
        dense_hpipm_solver_args->hpipm_args->warm_start = warmstart;

        // dense_hpipm_solver_args->hpipm_args->mu0 = 1e6;
        // dense_hpipm_solver_args->hpipm_args->res_g_max = 1e-6;  // < 1e-8 breaks hpipm
        // dense_hpipm_solver_args->hpipm_args->res_b_max = 1e-8;
        dense_hpipm_solver_args->hpipm_args->res_d_max = 1e-7;  // < 1e-7 breaks hpipm
        // dense_hpipm_solver_args->hpipm_args->res_m_max = 1e-8;
        // dense_hpipm_solver_args->hpipm_args->alpha_min = 1e-10;
        // dense_hpipm_solver_args->hpipm_args->cond_pred_corr = 1;

        break;
    case FULL_CONDENSING_QORE:
#ifdef ACADOS_WITH_QORE
        printf("\nFull condensing + QORE:\n\n");
        fcond_solver_args = (ocp_qp_full_condensing_solver_args *)args;
        dense_qp_qore_args *qore_solver_args = (dense_qp_qore_args *)fcond_solver_args->solver_args;

        qore_solver_args->max_iter = 1000;
        qore_solver_args->warm_start = warmstart;  // TODO(dimitris): ONLY WORKS COLD STARTED

        if (qore_solver_args->warm_start)
            qore_solver_args->warm_strategy = 0;

        break;
#endif
    case FULL_CONDENSING_QPOASES:
        printf("\nFull condensing + QPOASES:\n\n");
        fcond_solver_args = (ocp_qp_full_condensing_solver_args *)args;
        dense_qp_qpoases_args *qpoases_solver_args = (dense_qp_qpoases_args *)fcond_solver_args->solver_args;

        qpoases_solver_args->warm_start = warmstart;

        break;
    case PARTIAL_CONDENSING_OOQP:
        break;
    }

    ocp_qp_in *qp_in = create_ocp_qp_in(&dims);

    ocp_qp_out *qp_out = create_ocp_qp_out(&dims);

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

    // logs
    int nvars = 0;
    for (int ii = 0; ii < dims.N+1; ii++)
        nvars += dims.nx[ii] + dims.nu[ii];

    double *acado_sol;

    int *acado_iter_ptr;

    double *min_cpu_times = malloc(n_problems*sizeof(double));
    double *sol_error = malloc(n_problems*sizeof(double));
    int *iters = malloc(n_problems*sizeof(int));

    // info
    ocp_qp_info *info = (ocp_qp_info *)qp_out->misc;

    /************************************************
    * closed loop simulation
    ************************************************/

    ocp_qp_solver *qp_solver;

    for (int irun = 0; irun < n_rep; irun++)
    {
        qp_solver = ocp_qp_create(&plan, &dims, args);

        for (int indx = 0; indx < n_problems; indx++)
        {

            /************************************************
            * set up dynamics
            ************************************************/

            snprintf(str, sizeof(str), "Av%s_%d", suffix, indx);
            load_ptr(lib, str, (void **)&Av);

            snprintf(str, sizeof(str), "Bv_%d", indx);
            load_ptr(lib, str, (void **)&Bv);

            snprintf(str, sizeof(str), "b%s_%d", suffix, indx);
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

            snprintf(str, sizeof(str), "Qv%s_%d", suffix, indx);
            load_ptr(lib, str, (void **)&Qv);

            snprintf(str, sizeof(str), "Rv_%d", indx);
            load_ptr(lib, str, (void **)&Rv);

            snprintf(str, sizeof(str), "Sv%s_%d", suffix, indx);
            load_ptr(lib, str, (void **)&Sv);

            snprintf(str, sizeof(str), "q%s_%d", suffix, indx);
            load_ptr(lib, str, (void **)&q);

            snprintf(str, sizeof(str), "r%s_%d", suffix, indx);
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

            snprintf(str, sizeof(str), "lb%s_%d", suffix, indx);
            load_ptr(lib, str, (void **)&lb);

            snprintf(str, sizeof(str), "ub%s_%d", suffix, indx);
            load_ptr(lib, str, (void **)&ub);

            sum_nb = 0;

            for (int ii = 0; ii < N+1; ii++)
            {
                hlb[ii] = &lb[sum_nb];
                hub[ii] = &ub[sum_nb];
                sum_nb += dims.nb[ii];
            }

            /************************************************
            * get acado stats
            ************************************************/

            snprintf(str, sizeof(str), "acado_iter_%d", indx);
            load_ptr(lib, str, (void **)&acado_iter_ptr);

            snprintf(str, sizeof(str), "acado_sol%s_%d", suffix, indx);
            load_ptr(lib, str, (void **)&acado_sol);

            /************************************************
            * convert data and solve qp
            ************************************************/

            d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub,
                hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, qp_in);

            // print_ocp_qp_in(qp_in);

            int acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

            // print_ocp_qp_out(qp_out);

            sol_error[indx] = compare_with_acado_solution(N, nvars, qp_out, acado_sol);

            if (acados_return != ACADOS_SUCCESS)
            {
                if (acados_return == ACADOS_MINSTEP)
                    printf("QP SOLVER RETURNED MIN STEP STATUS (%d)\n", acados_return);
                else if (acados_return == ACADOS_MAXITER)
                    printf("QP SOLVER RETURNED MAX ITER STATUS (%d)\n", acados_return);
                else if (acados_return == ACADOS_FAILURE)
                    printf("QP SOLVER FAILED\n");
                else
                    printf("QP SOLVER RETURNED UNKNOWN FLAG\n");
                return -1;
            }

            /************************************************
            * print and save results
            ************************************************/

            printf("\n--> problem %d solved with acados in %d iterations and %f ms",
                indx, info->num_iter, info->total_time*1000);

            printf("\n--> problem %d solved with ACADO  in %d iterations and error is solution is %e\n\n",
                indx, *acado_iter_ptr, sol_error[indx]);

            if (irun == 0)
            {
                min_cpu_times[indx] = info->total_time;
                iters[indx] = info->num_iter;
            } else
            {
                if(iters[indx] != info->num_iter)
                {
                    printf("inconsistent number of iterations between runs\n");
                    exit(1);
                }
                if (min_cpu_times[indx] > info->total_time)
                    min_cpu_times[indx] = info->total_time;
            }

        }
        printf("\n------ end of run #%d---------\n", irun+1);
        free(qp_solver);
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

        char save_str[256];
        char *lib_str_no_ext = strndup(lib_str, strlen(lib_str) - strlen(".so"));

    if (auto_choose_acados_solver)
    {
        snprintf(save_str, sizeof(save_str), "%s_acados_cpu_times.txt", lib_str_no_ext);
        write_double_vector_to_txt(min_cpu_times, n_problems, save_str);
        snprintf(save_str, sizeof(save_str), "%s_acados_iters.txt", lib_str_no_ext);
        write_int_vector_to_txt(iters, n_problems, save_str);
        snprintf(save_str, sizeof(save_str), "%s_sol_error.txt", lib_str_no_ext);
        write_double_vector_to_txt(sol_error, n_problems, save_str);
    } else
    {
        #if 0  // to enable comparison with other solvers
        char custom_str[256] = "dense_hpipm";

        snprintf(save_str, sizeof(save_str), "%s_acados_%s_cpu_times.txt", lib_str_no_ext, custom_str);
        write_double_vector_to_txt(min_cpu_times, n_problems, save_str);
        snprintf(save_str, sizeof(save_str), "%s_acados_%s_iters.txt", lib_str_no_ext, custom_str);
        write_int_vector_to_txt(iters, n_problems, save_str);
        snprintf(save_str, sizeof(save_str), "%s_%s_sol_error.txt", lib_str_no_ext, custom_str);
        write_double_vector_to_txt(sol_error, n_problems, save_str);
        #endif
    }

    free(min_cpu_times);
    free(sol_error);
    free(iters);

    printf("\nacados runs with N2 = %d\n", N2);

    }  // end of loop over N_ins

    return 0;
}



void choose_solver(int N, char *lib_str, int *N2, int *warm_start, ocp_qp_solver_t *qp_solver)
{
    char *solver_str = strstr(lib_str, "solver");
    char *B_str_full = strstr(solver_str, "_B");
    char *warm_str_full = strstr(solver_str, "warmstart_");

    if (strstr(solver_str, "HPMPC") != NULL)
    {
        printf("HPMPC\t\tdetected\n");
        *qp_solver = PARTIAL_CONDENSING_HPMPC;
    }
    if (strstr(solver_str, "qpOASES") != NULL)
    {
        printf("qpOASES\t\tdetected\n");
        *qp_solver = FULL_CONDENSING_QPOASES;
    }
    if (strstr(solver_str, "qpDUNES") != NULL)
    {
        printf("qpDUNES\t\tdetected\n");
        *qp_solver = PARTIAL_CONDENSING_QPDUNES;
    }

    if (B_str_full != NULL)
    {
        int B_len = strlen(B_str_full) - strlen("_B_warmstart_X.so");
        char *B_str = strndup(B_str_full+2, B_len);
        int B = 0;
        for (int ii = 0; ii < B_len; ii++)
            B += (B_str[B_len-ii-1] - '0')*pow(10,ii);
        *N2 = B > 0 ? N/B : 0;
        printf("N2 = %d\t\tdetected\n", *N2);
    } else
    {
        *N2 = 0;
        printf("N2 \t\tnot detected (setting to 0)\n");
    }
    if (warm_str_full != NULL)
    {
        *warm_start = warm_str_full[10] - '0';
        printf("warmstart = %d\tdetected\n", *warm_start);
    } else
    {
        printf("warmstart \tnot detected (setting to 0)\n");
        *warm_start = 0;
    }
    // exit(1);
}
/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


// TODO(dimitris): VALGRIND!!! (WITHOUT QPDUNES)

// external
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// acados
#include <acados/utils/print.h>

#include <acados/ocp_qp/ocp_qp_full_condensing.h>

// c interface
#include <acados_c/ocp_qp_interface.h>

// mass spring helper functions
// hard constraints
ocp_qp_xcond_solver_dims *create_ocp_qp_dims_mass_spring(ocp_qp_xcond_solver_config *config, int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(ocp_qp_dims *dims);
// soft constraints
// TODO
ocp_qp_dims *create_ocp_qp_dims_mass_spring_soft_constr(int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring_soft_constr(ocp_qp_dims *dims);

// #ifndef ACADOS_WITH_QPDUNES
#define ELIMINATE_X0
// #endif

#define GENERAL_CONSTRAINT_AT_TERMINAL_STAGE

// #define SOFT_CONSTRAINTS

#define NREP 1

int main() {
    printf("\n");
    printf("\n");
    printf("\n");
    printf(" mass spring example: acados ocp_qp solvers\n");
    printf("\n");
    printf("\n");
    printf("\n");

    /************************************************
     * set up dimensions
     ************************************************/

    int nx_ = 8;   // number of states (it has to be even for the mass-spring system test problem)

    int nu_ = 3;   // number of inputs (controllers) (it has to be at least 1 and
                   // at most nx_/2 for the mass-spring system test problem)

    int N = 15;    // horizon length
    int nb_ = 11;  // number of box constrained inputs and states
    int ng_ = 0;   // 4;  // number of general constraints

    #ifdef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
    int num_of_stages_equal_to_zero = 4;  // number of states to be enforced to zero at last stage
    int ngN = num_of_stages_equal_to_zero;
    #else
    int ngN = 0;
    #endif

    /************************************************
     * ocp qp solvers
     ************************************************/

    // choose values for N2 in partial condensing solvers
    int num_N2_values = 3;
    int N2_values[3] = {15,10,5};

    int ii_max = 8;

    #ifndef ACADOS_WITH_HPMPC
    ii_max--;
    #endif
    #ifndef ACADOS_WITH_QPDUNES
    ii_max--;
    #endif
    #ifndef ACADOS_WITH_QORE
    ii_max--;
    #endif
    #ifndef ACADOS_WITH_QPOASES
    ii_max--;
    #endif
    #ifndef ACADOS_WITH_OOQP
    ii_max--;
    ii_max--;
    #endif
    #ifndef ACADOS_WITH_OSQP
    ii_max--;
    #endif

    // choose ocp qp solvers
    ocp_qp_solver_t ocp_qp_solvers[] =
    {
        PARTIAL_CONDENSING_HPIPM,

        FULL_CONDENSING_HPIPM,

        #ifdef ACADOS_WITH_HPMPC
        PARTIAL_CONDENSING_HPMPC,
        #endif

        // #ifdef ACADOS_WITH_QPDUNES
        // PARTIAL_CONDENSING_QPDUNES,
        // #endif

        #ifdef ACADOS_WITH_QORE
        FULL_CONDENSING_QORE,
        #endif

        #ifdef ACADOS_WITH_QPOASES
        FULL_CONDENSING_QPOASES,
        #endif

        #ifdef ACADOS_WITH_OOQP
        PARTIAL_CONDENSING_OOQP,
        FULL_CONDENSING_OOQP,
        #endif

        #ifdef ACADOS_WITH_OSQP
        PARTIAL_CONDENSING_OSQP,
        #endif
    };


    /************************************************
     * simulations
     ************************************************/

    for (int ii = 0; ii < ii_max; ii++)
    {
        // plan
        ocp_qp_solver_plan_t plan;
        plan.qp_solver = ocp_qp_solvers[ii];

        // config
        ocp_qp_xcond_solver_config *config = ocp_qp_xcond_solver_config_create(plan);

        // dims
#ifdef SOFT_CONSTRAINTS
        // TODO
        ocp_qp_dims *qp_dims = create_ocp_qp_dims_mass_spring_soft_constr(N, nx_, nu_, nb_, ng_, ngN);
#else
        ocp_qp_xcond_solver_dims *qp_dims = create_ocp_qp_dims_mass_spring(config, N, nx_, nu_, nb_, ng_, ngN);
#endif

        // qp in
#ifdef SOFT_CONSTRAINTS
        // TODO
        ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring_soft_constr(qp_dims->orig_dims);
#else
        ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring(qp_dims->orig_dims);
#endif

        // qp out
        ocp_qp_out *qp_out = ocp_qp_out_create(qp_dims->orig_dims);


        ocp_qp_xcond_solver_opts *opts = ocp_qp_xcond_solver_opts_create(config, qp_dims);
//        bool ok = false;
//        if (ok == true) config++; // dummy command to shut up Werror in Release

#ifdef ACADOS_WITH_QPDUNES
        int clipping;
#endif

#ifdef SOFT_CONSTRAINTS
        double mu0 = 1e2;
#endif

        for (int jj = 0; jj < num_N2_values; jj++)
        {
            int N2, max_iter, warm_start;

            // NOTE(nielsvd): needs to be implemented using the acados_c/options.h interface
            switch (plan.qp_solver)
            {
                case PARTIAL_CONDENSING_HPIPM:
                    printf("\nPartial condensing + HPIPM (N2 = %d):\n\n", N2);

                    N2 = N2_values[jj];
                    config->opts_set(config, opts, "cond_N", &N2);
                    max_iter = 30;
                    config->opts_set(config, opts, "iter_max", &max_iter);
//                    ok = set_option_double(opts, "sparse_hpipm.res_g_max", 1e-8);
//                    assert(ok = true && "specified option not found!");
//                    ok = set_option_double(opts, "sparse_hpipm.res_b_max", 1e-8);
//                    assert(ok = true && "specified option not found!");
//                    ok = set_option_double(opts, "sparse_hpipm.res_d_max", 1e-8);
//                    assert(ok = true && "specified option not found!");
//                    ok = set_option_double(opts, "sparse_hpipm.res_m_max", 1e-8);
//                    assert(ok = true && "specified option not found!");
#ifdef SOFT_CONSTRAINTS
                    mu0 = 1e2;
                    config->opts_set(config, opts, "mu0", &mu0);
#endif

                    break;
#ifdef ACADOS_WITH_HPMPC
                case PARTIAL_CONDENSING_HPMPC:
                    printf("\nPartial condensing + HPMPC (N2 = %d):\n\n", N2);

                    N2 = N2_values[jj];
                    config->opts_set(config, opts, "cond_N", &N2);
                    max_iter = 30;
                    config->opts_set(config, opts, "iter_max", &max_iter);
                    break;
#endif
// #ifdef ACADOS_WITH_QPDUNES
//                 case PARTIAL_CONDENSING_QPDUNES:
//                     printf("\nPartial condensing + qpDUNES (N2 = %d):\n\n", N2);
//                 #ifdef ELIMINATE_X0
//                     assert(1==0 && "qpDUNES does not support ELIMINATE_X0 flag!");
//                 #endif

//                 #ifdef GENERAL_CONSTRAINT_AT_TERMINAL_STAGE
//                     clipping = 0;
//                     config->opts_set(config, opts, "clipping", &clipping);
//                 #else
//                     if (N2 == N)
//                     {
//                         clipping = 1;
//                         config->opts_set(config, opts, "clipping", &clipping);
//                     } else
//                     {
//                         clipping = 0;
//                         config->opts_set(config, opts, "clipping", &clipping);
//                     }
//                 #endif
//                     warm_start = 0;
//                     config->opts_set(config, opts, "warm_start", &warm_start);

//                     N2 = N2_values[jj];
//                     config->opts_set(config, opts, "cond_N", &N2);
//                     break;
// #endif
                case FULL_CONDENSING_HPIPM:
                    printf("\nFull condensing + HPIPM:\n\n");
                    // default options
                    break;
#ifdef ACADOS_WITH_QORE
                case FULL_CONDENSING_QORE:
                    printf("\nFull condensing + QORE:\n\n");
                    // default options
                    break;
#endif
#ifdef ACADOS_WITH_QPOASES
                case FULL_CONDENSING_QPOASES:
                    printf("\nFull condensing + QPOASES:\n\n");
                    warm_start = 0;
                    config->opts_set(config, opts, "warm_start", &warm_start);
                    break;
#endif
#ifdef ACADOS_WITH_OOQP
                case PARTIAL_CONDENSING_OOQP:
                    printf("\nPartial condensing + OOQP (N2 = %d):\n\n", N2);
                    N2 = N2_values[jj];
                    config->opts_set(config, opts, "cond_N", &N2);
                    break;

                case FULL_CONDENSING_OOQP:
                    printf("\nFull condensing + OOQP:\n\n");
                    break;
#endif
#ifdef ACADOS_WITH_OSQP
                case PARTIAL_CONDENSING_OSQP:
                    printf("\nPartial condensing + OSQP (N2 = %d):\n\n", N2);
                    N2 = N2_values[jj];
                    config->opts_set(config, opts, "cond_N", &N2);
                    break;

#endif
                case INVALID_QP_SOLVER:
                    printf("\nInvalid QP solver\n\n");

            }

//            ocp_qp_full_condensing_dims *xcond_dims = qp_dims->xcond_dims;
//            printf("%d %d\n", xcond_dims->fcond_dims->nv, xcond_dims->fcond_dims->ne);
//            printf("\nbefore solver create\n");
            ocp_qp_solver *qp_solver = ocp_qp_create(config, qp_dims, opts);
//            printf("\nafter solver create\n");
//            exit(1);

            int acados_return = 0;

            qp_info *info = (qp_info *) qp_out->misc;
            qp_info min_info;

            // print_ocp_qp_in(qp_in);

            // run QP solver NREP times and record min timings
            for (int rep = 0; rep < NREP; rep++)
            {
                acados_return += ocp_qp_solve(qp_solver, qp_in, qp_out);

                if (rep == 0)
                {
                    min_info.num_iter = info->num_iter;
                    min_info.total_time = info->total_time;
                    min_info.condensing_time = info->condensing_time;
                    min_info.solve_QP_time = info->solve_QP_time;
                    min_info.interface_time = info->interface_time;
                }
                else
                {
                    assert(min_info.num_iter == info->num_iter && "QP solver not cold started!");

                    if (info->total_time < min_info.total_time)
                        min_info.total_time = info->total_time;
                    if (info->condensing_time < min_info.condensing_time)
                        min_info.condensing_time = info->condensing_time;
                    if (info->solve_QP_time < min_info.solve_QP_time)
                        min_info.solve_QP_time = info->solve_QP_time;
                    if (info->interface_time < min_info.interface_time)
                        min_info.interface_time = info->interface_time;
                }
            }

            /************************************************
             * compute infinity norm of residuals
             ************************************************/

            double res[4];
            ocp_qp_inf_norm_residuals(qp_dims->orig_dims, qp_in, qp_out, res);

            double max_res = 0.0;
            for (int ii = 0; ii < 4; ii++)
                max_res = (res[ii] > max_res) ? res[ii] : max_res;

            /************************************************
             * print solutions and stats
             ************************************************/

            print_ocp_qp_out(qp_out);

            printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);

            print_qp_info(&min_info);

            free(qp_solver);

            // NOTE(dimitris): run dense solver once (sparse solvers multiple times for N2 values)
            if (plan.qp_solver >= FULL_CONDENSING_HPIPM) break;
        }
        ocp_qp_xcond_solver_config_free(config);
        ocp_qp_xcond_solver_opts_free(opts);
        ocp_qp_in_free(qp_in);
        ocp_qp_out_free(qp_out);
    }

    printf("\nsuccess!\n\n");
}

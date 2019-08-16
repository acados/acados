/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


#include <stdlib.h>

#include "acados_c/options_interface.h"

#include "acados/dense_qp/dense_qp_hpipm.h"
#ifdef ACADOS_WITH_QORE
#include "acados/dense_qp/dense_qp_qore.h"
#endif
#ifdef ACADOS_WITH_QPOASES
#include "acados/dense_qp/dense_qp_qpoases.h"
#endif
#include "acados/ocp_qp/ocp_qp_full_condensing_solver.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#ifdef ACADOS_WITH_HPMPC
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#endif
#ifdef ACADOS_WITH_OOQP
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "acados/dense_qp/dense_qp_ooqp.h"
#endif
#ifdef ACADOS_WITH_OSQP
#include "acados/ocp_qp/ocp_qp_osqp.h"
#endif
#ifdef ACADOS_WITH_QPDUNES
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#endif
#include "acados/ocp_qp/ocp_qp_partial_condensing_solver.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_nlp/ocp_nlp_sqp_rti.h"
#include "acados/utils/strsep.h"

static bool is_qp_solver(const char *qp_solver_name)
{
    return !strcmp(qp_solver_name, "hpipm")
        || !strcmp(qp_solver_name, "sparse_hpipm")
        || !strcmp(qp_solver_name, "condensing_hpipm")
        || !strcmp(qp_solver_name, "hpmpc")
        || !strcmp(qp_solver_name, "ooqp")
        || !strcmp(qp_solver_name, "qpdunes")
        || !strcmp(qp_solver_name, "qpoases")
        || !strcmp(qp_solver_name, "qore");
}

bool set_option_int(void *args_, const char *option, const int value)
{
    char *token;
    char option_cpy[MAX_STR_LEN];
    strcpy(option_cpy, option);
    char *ptr_to_option_cpy = &option_cpy[0];
    token = strsep_acados(&ptr_to_option_cpy, ".");
    while (token)
    {
        // Linear search since the number of options is small.
        if (!strcmp(token, "sqp"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_nlp_sqp_opts *args = (ocp_nlp_sqp_opts *) args_;
            if (!strcmp(token, "max_iter")) {
                args->max_iter = value;
            }
            else if (is_qp_solver(token))
            {
                token[strlen(token)] = '.';  // this effectively concatenates the strings again
                return set_option_int(args->qp_solver_opts, token, value);
            }
        }
        else if (!strcmp(token, "rti"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_nlp_sqp_rti_opts *args = (ocp_nlp_sqp_rti_opts *) args_;
            if (is_qp_solver(token))
            {
                token[strlen(token)] = '.';  // this effectively concatenates the strings again
                return set_option_int(args->qp_solver_opts, token, value);
            }
        }
        else if (!strcmp(token, "sparse_hpipm") || !strcmp(token, "hpipm"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_hpipm_opts *args = (ocp_qp_hpipm_opts *) sparse_args->qp_solver_opts;
            ocp_qp_partial_condensing_opts *pcond_opts =
                (ocp_qp_partial_condensing_opts *) sparse_args->pcond_opts;
            if (!strcmp(token, "max_iter"))
                args->hpipm_opts->iter_max = value;
            else if (!strcmp(token, "max_stat"))
                args->hpipm_opts->stat_max = value;
            else if (!strcmp(token, "cond_N"))
                pcond_opts->N2 = value;
            else
                return false;
        }
        else if (!strcmp(token, "condensing_hpipm"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            dense_qp_hpipm_opts *args = (dense_qp_hpipm_opts *) sparse_args->qp_solver_opts;
            if (!strcmp(token, "max_iter"))
                args->hpipm_opts->iter_max = value;
            else if (!strcmp(token, "max_stat"))
                args->hpipm_opts->stat_max = value;
            else
                return false;
#ifdef ACADOS_WITH_HPMPC
        }
        else if (!strcmp(token, "hpmpc"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_hpmpc_opts *args = (ocp_qp_hpmpc_opts *) sparse_args->qp_solver_opts;
            ocp_qp_partial_condensing_opts *pcond_opts =
                (ocp_qp_partial_condensing_opts *) sparse_args->pcond_opts;
            if (!strcmp(token, "max_iter"))
                args->max_iter = value;
            else if (!strcmp(token, "warm_start"))
                args->warm_start = value;
            // NOTE(dimitris): HPMPC partial condesing has a bug, using hpipm partial condensing
            // instead
            else if (!strcmp(token, "cond_N"))
                pcond_opts->N2 = value;
            // partial tightening
            else if (!strcmp(token, "N"))
                args->N = value;
            else if (!strcmp(token, "M"))
                args->M = value;
            else
                return false;
#endif
#ifdef ACADOS_WITH_OOQP
        }
        else if (!strcmp(token, "sparse_ooqp"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_ooqp_opts *args = (ocp_qp_ooqp_opts *) sparse_args->qp_solver_opts;
            ocp_qp_partial_condensing_opts *pcond_opts =
                (ocp_qp_partial_condensing_opts *) sparse_args->pcond_opts;
            if (!strcmp(token, "print_level"))
                args->printLevel = value;
            else if (!strcmp(token, "cond_N"))
                pcond_opts->N2 = value;
            else
                return false;
        }
        else if (!strcmp(token, "condensing_ooqp"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            dense_qp_ooqp_opts *args = (dense_qp_ooqp_opts *) sparse_args->qp_solver_opts;
            if (!strcmp(token, "print_level"))
                args->printLevel = value;
            else
                return false;
#endif
#ifdef ACADOS_WITH_OSQP
        }
        else if (!strcmp(token, "sparse_osqp"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_osqp_opts *args = (ocp_qp_osqp_opts *) sparse_args->qp_solver_opts;
            ocp_qp_partial_condensing_opts *pcond_opts =
                (ocp_qp_partial_condensing_opts *) sparse_args->pcond_opts;
            if (!strcmp(token, "cond_N"))
                pcond_opts->N2 = value;
            else
                return false;
#endif
#ifdef ACADOS_WITH_QPDUNES
        }
        else if (!strcmp(token, "qpdunes"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_qpdunes_opts *args = (ocp_qp_qpdunes_opts *) sparse_args->qp_solver_opts;
            ocp_qp_partial_condensing_opts *pcond_opts =
                (ocp_qp_partial_condensing_opts *) sparse_args->pcond_opts;
            if (!strcmp(token, "print_level"))
            {
                args->options.printLevel = value;
            }
            else if (!strcmp(token, "warm_start"))
            {
                args->warmstart = value;
            }
            else if (!strcmp(token, "max_iter"))
            {
                args->options.maxIter = value;
            }
            else if (!strcmp(token, "cond_N"))
            {
                pcond_opts->N2 = value;
            }
            else if (!strcmp(token, "clipping"))
            {
                if (value == 1)
                {
                    args->stageQpSolver = QPDUNES_WITH_CLIPPING;
                    args->options.lsType = QPDUNES_LS_ACCELERATED_GRADIENT_BISECTION_LS;
                }
                else
                {
                    args->stageQpSolver = QPDUNES_WITH_QPOASES;
                    args->options.lsType = QPDUNES_LS_HOMOTOPY_GRID_SEARCH;
                }
            }
            else
            {
                return false;
            }
#endif
#ifdef ACADOS_WITH_QORE
        }
        else if (!strcmp(token, "qore"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_full_condensing_solver_opts *cond_opts =
                (ocp_qp_full_condensing_solver_opts *) args_;
            dense_qp_qore_opts *args = (dense_qp_qore_opts *) cond_opts->qp_solver_opts;
            if (!strcmp(token, "print_freq"))
                args->print_freq = value;
            else if (!strcmp(token, "warm_start"))
                args->warm_start = value;
            else if (!strcmp(token, "warm_strategy"))
                args->warm_strategy = value;
            else if (!strcmp(token, "nsmax"))
                args->nsmax = value;
            else if (!strcmp(token, "hot_start"))
                args->hot_start = value;
            else if (!strcmp(token, "max_iter"))
                args->max_iter = value;
            else
                return false;
#endif
#ifdef ACADOS_WITH_QPOASES
        }
        else if (!strcmp(token, "qpoases"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_full_condensing_solver_opts *cond_opts =
                (ocp_qp_full_condensing_solver_opts *) args_;
            dense_qp_qpoases_opts *args = (dense_qp_qpoases_opts *) cond_opts->qp_solver_opts;
            if (!strcmp(token, "max_iter"))
                args->max_nwsr = value;
            else if (!strcmp(token, "warm_start"))
                args->warm_start = value;
            else
                return false;
#endif
        }
        else
        {
            return false;
        }
        token = strsep_acados(&ptr_to_option_cpy, ".");
    }
    return true;
}


bool set_option_double(void *args_, const char *option, const double value)
{
    char *token;
    char option_cpy[MAX_STR_LEN];
    strcpy(option_cpy, option);
    char *ptr_to_option_cpy = &option_cpy[0];
    while ((token = strsep_acados(&ptr_to_option_cpy, ".")))
    {
        // Linear search since the number of options is small.
        if (!strcmp(token, "sparse_hpipm"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_hpipm_opts *args = (ocp_qp_hpipm_opts *) sparse_args->qp_solver_opts;
            if (!strcmp(token, "res_g_max"))
                args->hpipm_opts->res_g_max = value;
            else if (!strcmp(token, "res_b_max"))
                args->hpipm_opts->res_b_max = value;
            else if (!strcmp(token, "res_d_max"))
                args->hpipm_opts->res_d_max = value;
            else if (!strcmp(token, "res_m_max"))
                args->hpipm_opts->res_m_max = value;
            else if (!strcmp(token, "alpha_min"))
                args->hpipm_opts->alpha_min = value;
            else if (!strcmp(token, "mu0"))
                args->hpipm_opts->mu0 = value;
            else
                return false;
        }
        else if (!strcmp(token, "condensing_hpipm"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            dense_qp_hpipm_opts *args = (dense_qp_hpipm_opts *) sparse_args->qp_solver_opts;
            if (!strcmp(token, "res_g_max"))
                args->hpipm_opts->res_g_max = value;
            else if (!strcmp(token, "res_b_max"))
                args->hpipm_opts->res_b_max = value;
            else if (!strcmp(token, "res_d_max"))
                args->hpipm_opts->res_d_max = value;
            else if (!strcmp(token, "res_m_max"))
                args->hpipm_opts->res_m_max = value;
            else if (!strcmp(token, "alpha_min"))
                args->hpipm_opts->alpha_min = value;
            else if (!strcmp(token, "mu0"))
                args->hpipm_opts->mu0 = value;
            else
                return false;
#ifdef ACADOS_WITH_HPMPC
        }
        else if (!strcmp(token, "hpmpc"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_hpmpc_opts *args = (ocp_qp_hpmpc_opts *) sparse_args->qp_solver_opts;
            if (!strcmp(token, "tol"))
                args->tol = value;
            else if (!strcmp(token, "mu0"))
                args->mu0 = value;
            // partial tightening
            else if (!strcmp(token, "sigma_mu"))
                args->sigma_mu = value;
            else
                return false;
#endif
#ifdef ACADOS_WITH_QPDUNES
        }
        else if (!strcmp(token, "qpdunes"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_partial_condensing_solver_opts *sparse_args =
                (ocp_qp_partial_condensing_solver_opts *) args_;
            ocp_qp_qpdunes_opts *args = (ocp_qp_qpdunes_opts *) sparse_args->qp_solver_opts;
            if (!strcmp(token, "tolerance"))
                args->options.stationarityTolerance = value;
            else
                return false;
#endif
#ifdef ACADOS_WITH_QPOASES
        }
        else if (!strcmp(token, "qpoases"))
        {
            token = strsep_acados(&ptr_to_option_cpy, ".");
            ocp_qp_full_condensing_solver_opts *cond_opts =
                (ocp_qp_full_condensing_solver_opts *) args_;
            dense_qp_qpoases_opts *args = (dense_qp_qpoases_opts *) cond_opts->qp_solver_opts;
            if (!strcmp(option, "max_cputime"))
                args->max_cputime = value;
            else
                return false;
#endif
        }
        else
        {
            return false;
        }
        token = strsep_acados(&ptr_to_option_cpy, ".");
    }
    return true;
}

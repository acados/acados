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


#include "acados/ocp_nlp/ocp_nlp_sqp.h"

// external
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if defined(ACADOS_WITH_OPENMP)
#include <omp.h>
#endif

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados_c/ocp_qp_interface.h"



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_sqp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_sqp_opts);

    size += ocp_nlp_opts_calculate_size(config, dims);

    return size;
}



void *ocp_nlp_sqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_opts);

    opts->nlp_opts = ocp_nlp_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_sqp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_sqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    // int ii;

    // this first !!!
    ocp_nlp_opts_initialize_default(config, dims, nlp_opts);

    // SQP opts
    opts->max_iter = 20;
    opts->tol_stat = 1e-8;
    opts->tol_eq   = 1e-8;
    opts->tol_ineq = 1e-8;
    opts->tol_comp = 1e-8;
    opts->tol_unbounded = -1e10;
    opts->tol_min_step_norm = 1e-12;

    opts->ext_qp_res = 0;

    opts->qp_warm_start = 0;
    opts->warm_start_first_qp = false;
    opts->rti_phase = 0;
    opts->eval_residual_at_max_iter = false;

    // funnel method opts
    opts->funnel_initialization_increase_factor = 15.0;
    opts->funnel_initialization_upper_bound = 1.0;
    opts->funnel_sufficient_decrease_factor = 0.9;
    opts->funnel_kappa = 0.9;
    opts->funnel_fraction_switching_condition = 1e-3;
    opts->funnel_initial_penalty_parameter = 1.0;
    opts->funnel_penalty_contraction = 0.5;
    opts->funnel_penalty_eta = 1e-6;
    opts->funnel_type_switching_condition = false; // use ipopt/gould type of switching

    // overwrite default submodules opts
    // qp tolerance
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", &opts->tol_stat);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", &opts->tol_eq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", &opts->tol_ineq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", &opts->tol_comp);

    return;
}



void ocp_nlp_sqp_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_update(config, dims, nlp_opts);

    return;
}



void ocp_nlp_sqp_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    int ii;

    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (ii=0; ii<module_length; ii++)
            module[ii] = field[ii];
        module[module_length] = '\0'; // add end of string
        ptr_module = module;
    }

    // pass options to QP module
    if ( ptr_module!=NULL && (!strcmp(ptr_module, "qp")) )
    {
        ocp_nlp_opts_set(config, nlp_opts, field, value);

        if (!strcmp(field, "qp_warm_start"))
        {
            int* i_ptr = (int *) value;
            opts->qp_warm_start = *i_ptr;
        }
    }
    else // nlp opts
    {
        if (!strcmp(field, "max_iter"))
        {
            int* max_iter = (int *) value;
            opts->max_iter = *max_iter;
        }
        else if (!strcmp(field, "tol_stat"))
        {
            double* tol_stat = (double *) value;
            opts->tol_stat = *tol_stat;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", value);
        }
        else if (!strcmp(field, "tol_eq"))
        {
            double* tol_eq = (double *) value;
            opts->tol_eq = *tol_eq;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", value);
        }
        else if (!strcmp(field, "tol_ineq"))
        {
            double* tol_ineq = (double *) value;
            opts->tol_ineq = *tol_ineq;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", value);
        }
        else if (!strcmp(field, "tol_comp"))
        {
            double* tol_comp = (double *) value;
            opts->tol_comp = *tol_comp;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", value);
        }
        else if (!strcmp(field, "tol_min_step_norm"))
        {
            double* tol_min_step_norm = (double *) value;
            opts->tol_min_step_norm = *tol_min_step_norm;
        }
        else if (!strcmp(field, "ext_qp_res"))
        {
            int* ext_qp_res = (int *) value;
            opts->ext_qp_res = *ext_qp_res;
        }
        else if (!strcmp(field, "warm_start_first_qp"))
        {
            bool* warm_start_first_qp = (bool *) value;
            opts->warm_start_first_qp = *warm_start_first_qp;
        }
        else if (!strcmp(field, "rti_phase"))
        {
            int* rti_phase = (int *) value;
            if (*rti_phase < 0 || *rti_phase > 0) {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for rti_phase field.");
                printf("possible values are: 0\n");
                exit(1);
            }
            opts->rti_phase = *rti_phase;
        }
        else if (!strcmp(field, "eval_residual_at_max_iter"))
        {
            bool* eval_residual_at_max_iter = (bool *) value;
            opts->eval_residual_at_max_iter = *eval_residual_at_max_iter;
        }
        else if (!strcmp(field, "funnel_initialization_increase_factor"))
        {
            double* funnel_initialization_increase_factor = (double *) value;
            if (*funnel_initialization_increase_factor <= 1.0)
            {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for funnel_initialization_increase_factor field, need double > 1, got %f.", *funnel_initialization_increase_factor);
                exit(1);
            }
            opts->funnel_initialization_increase_factor = *funnel_initialization_increase_factor;
        }
        else if (!strcmp(field, "funnel_initialization_upper_bound"))
        {
            double* funnel_initialization_upper_bound = (double *) value;
            if (*funnel_initialization_upper_bound <= 0.0)
            {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for funnel_initialization_upper_bound field, need double > 0, got %f.", *funnel_initialization_upper_bound);
                exit(1);
            }
            opts->funnel_initialization_upper_bound = *funnel_initialization_upper_bound;
        }
        else if (!strcmp(field, "funnel_sufficient_decrease_factor"))
        {
            double* funnel_sufficient_decrease_factor = (double *) value;
            if (*funnel_sufficient_decrease_factor <= 0.0 || *funnel_sufficient_decrease_factor >= 1.0)
            {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for funnel_sufficient_decrease_factor field, need double in (0,1), got %f.", *funnel_sufficient_decrease_factor);
                exit(1);
            }
            opts->funnel_sufficient_decrease_factor = *funnel_sufficient_decrease_factor;
        }
        else if (!strcmp(field, "funnel_kappa"))
        {
            double* funnel_kappa = (double *) value;
            if (*funnel_kappa <= 0.0 || *funnel_kappa >= 1.0)
            {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for funnel_kappa field, need double in (0,1), got %f.", *funnel_kappa);
                exit(1);
            }
            opts->funnel_kappa = *funnel_kappa;
        }
        else if (!strcmp(field, "funnel_fraction_switching_condition"))
        {
            double* funnel_fraction_switching_condition = (double *) value;
            if (*funnel_fraction_switching_condition <= 0.0 || *funnel_fraction_switching_condition >= 1.0)
            {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for funnel_fraction_switching_condition field, need double in (0,1), got %f.", *funnel_fraction_switching_condition);
                exit(1);
            }
            opts->funnel_fraction_switching_condition = *funnel_fraction_switching_condition;
        }
        else if (!strcmp(field, "funnel_initial_penalty_parameter"))
        {
            double* funnel_initial_penalty_parameter = (double *) value;
            if (*funnel_initial_penalty_parameter < 0.0 || *funnel_initial_penalty_parameter > 1.0)
            {
                printf("\nerror: ocp_nlp_sqp_opts_set: invalid value for funnel_initial_penalty_parameter field, need double in [0,1], got %f.", *funnel_initial_penalty_parameter);
                exit(1);
            }
            opts->funnel_initial_penalty_parameter = *funnel_initial_penalty_parameter;
        }
        else
        {
            ocp_nlp_opts_set(config, nlp_opts, field, value);
        }
    }

    return;

}



void ocp_nlp_sqp_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_set_at_stage(config, nlp_opts, stage, field, value);

    return;

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_sqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_sqp_memory);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // primal step norm
    if (opts->nlp_opts->log_primal_step_norm)
    {
        size += opts->max_iter*sizeof(double);
    }
    // stat
    int stat_m = opts->max_iter+1;
    int stat_n = 7;
    if (opts->ext_qp_res)
        stat_n += 4;
    size += stat_n*stat_m*sizeof(double);

    size += 3*8;  // align

    make_int_multiple_of(8, &size);

    return size;
}

void *ocp_nlp_sqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_dynamics_config **dynamics = config->dynamics;
    // ocp_nlp_cost_config **cost = config->cost;
    // ocp_nlp_constraints_config **constraints = config->constraints;

    char *c_ptr = (char *) raw_memory;

    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_sqp_memory *mem = (ocp_nlp_sqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_memory);

    align_char_to(8, &c_ptr);

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, nlp_opts, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // primal step norm
    if (opts->nlp_opts->log_primal_step_norm)
    {
        mem->primal_step_norm = (double *) c_ptr;
        c_ptr += opts->max_iter*sizeof(double);
    }

    // stat
    mem->stat = (double *) c_ptr;
    mem->stat_m = opts->max_iter+1;
    mem->stat_n = 7;
    if (opts->ext_qp_res)
        mem->stat_n += 4;
    c_ptr += mem->stat_m*mem->stat_n*sizeof(double);

    mem->status = ACADOS_READY;

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_sqp_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}

/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_sqp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    // sqp
    size += sizeof(ocp_nlp_sqp_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    if (opts->ext_qp_res)
    {
        // qp res
        size += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        size += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    return size;
}



static void ocp_nlp_sqp_cast_workspace(ocp_nlp_config *config, ocp_nlp_dims *dims,
         ocp_nlp_sqp_opts *opts, ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_workspace *work)
{
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_workspace);

    // nlp
    work->nlp_work = ocp_nlp_workspace_assign(config, dims, nlp_opts, nlp_mem, c_ptr);
    c_ptr += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    if (opts->ext_qp_res)
    {
        // qp res
        work->qp_res = ocp_qp_res_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        work->qp_res_ws = ocp_qp_res_workspace_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    assert((char *) work + ocp_nlp_sqp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}


/* Helper functions */


#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
static void ocp_nlp_sqp_dump_qp_in_to_file(ocp_qp_in *qp_in, int sqp_iter, int soc)
{
    char filename[100];
    if (soc > 0)
        sprintf(filename, "soc_qp_in_%d.txt", sqp_iter);
    else
        sprintf(filename, "qp_in_%d.txt", sqp_iter);
    FILE *out_file = fopen(filename, "w");
    print_ocp_qp_in_to_file(out_file, qp_in);
    fclose(out_file);
}


static void ocp_nlp_sqp_dump_qp_out_to_file(ocp_qp_out *qp_out, int sqp_iter, int soc)
{
    char filename[100];
    if (soc > 0)
        sprintf(filename, "soc_qp_out_%d.txt", sqp_iter);
    else
        sprintf(filename, "qp_out_%d.txt", sqp_iter);
    FILE *out_file = fopen(filename, "w");
    print_ocp_qp_out_to_file(out_file, qp_out);
    fclose(out_file);
}
#endif


static bool ocp_nlp_soc_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
            ocp_nlp_out *nlp_out, ocp_nlp_sqp_opts *opts, ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_workspace *work, int sqp_iter)
{
    int ii;
    int N = dims->N;

    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_qp_in *qp_in = nlp_mem->qp_in;
    ocp_qp_out *qp_out = nlp_mem->qp_out;
    qp_info *qp_info_;
    // NOTE: following Waechter2006:
    // Do SOC
    // 1. if "the first trial step size alpha_k,0 has been rejected and
    // 2. if the infeasibility would have increased when accepting the previous step
    // NOTE: the "and" is interpreted as an "or" in the current implementation

    // preliminary line search
    mem->alpha = ocp_nlp_line_search(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, 1, sqp_iter);
    if (mem->alpha >= 1.0)
    {
        return false; // do_line_search;
    }

    // Second Order Correction (SOC): following Nocedal2006: p.557, eq. (18.51) -- (18.56)
    // Paragraph: APPROACH III: S l1 QP (SEQUENTIAL l1 QUADRATIC PROGRAMMING),
    // Section 18.8 TRUST-REGION SQP METHODS
    //   - just no trust region radius here.
    if (nlp_opts->print_level > 0)
        printf("ocp_nlp_sqp: performing SOC, since alpha %e in prelim. line search\n\n", mem->alpha);
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;
    // int *nv = dims->nv;
    // int *ni = dims->ni;

    /* evaluate constraints & dynamics at new step */
    // NOTE: setting up the new iterate and evaluating is not needed here,
    //   since this evaluation was perfomed just before this call in the early terminated line search.

    // NOTE: similar to ocp_nlp_evaluate_merit_fun
    // update QP rhs
    // d_i = c_i(x_k + p_k) - \nabla c_i(x_k)^T * p_k
    struct blasfeo_dvec *tmp_fun_vec;

    for (ii = 0; ii <= N; ii++)
    {
        if (ii < N)
        {
            // b -- dynamics
            tmp_fun_vec = config->dynamics[ii]->memory_get_fun_ptr(nlp_mem->dynamics[ii]);
            // add - \nabla c_i(x_k)^T * p_k
            // c_i = f(x_k, u_k) - x_{k+1} (see dynamics module)
            blasfeo_dgemv_t(nx[ii]+nu[ii], nx[ii+1], -1.0, qp_in->BAbt+ii, 0, 0,
                            qp_out->ux+ii, 0, -1.0, tmp_fun_vec, 0, qp_in->b+ii, 0);
            // NOTE: not sure why it is - tmp_fun_vec here!
            blasfeo_dvecad(nx[ii+1], 1.0, qp_out->ux+ii+1, nu[ii+1], qp_in->b+ii, 0);
        }

        /* INEQUALITIES */
        // d -- constraints
        tmp_fun_vec = config->constraints[ii]->memory_get_fun_ptr(nlp_mem->constraints[ii]);
        /* SOC for bounds can be skipped (because linear) */
        // NOTE: SOC can also be skipped for truely linear constraint, i.e. ng of nlp,
        //      now using ng of QP = (nh+ng)

        // upper & lower
        blasfeo_dveccp(ng[ii], tmp_fun_vec, nb[ii], qp_in->d+ii, nb[ii]); // lg
        blasfeo_dveccp(ng[ii], tmp_fun_vec, 2*nb[ii]+ng[ii], qp_in->d+ii, 2*nb[ii]+ng[ii]); // ug
        // general linear / linearized!
        // tmp_ni = D * u + C * x
        blasfeo_dgemv_t(nu[ii]+nx[ii], ng[ii], 1.0, qp_in->DCt+ii, 0, 0, qp_out->ux+ii, 0,
                        0.0, &work->nlp_work->tmp_ni, 0, &work->nlp_work->tmp_ni, 0);
        // d[nb:nb+ng] += tmp_ni (lower)
        blasfeo_dvecad(ng[ii], 1.0, &work->nlp_work->tmp_ni, 0, qp_in->d+ii, nb[ii]);
        // d[nb:nb+ng] -= tmp_ni
        blasfeo_dvecad(ng[ii], -1.0, &work->nlp_work->tmp_ni, 0, qp_in->d+ii, 2*nb[ii]+ng[ii]);

        // add slack contributions
        // d[nb:nb+ng] += slack[idx]
        // qp_in->idxs_rev
        for (int j = 0; j < nb[ii]+ng[ii]; j++)
        {
            int slack_index = qp_in->idxs_rev[ii][j];
            if (slack_index >= 0)
            {
                // add slack contribution for lower and upper constraint
                // lower
                BLASFEO_DVECEL(qp_in->d+ii, j) -=
                        BLASFEO_DVECEL(qp_out->ux+ii, slack_index+nx[ii]+nu[ii]);
                // upper
                BLASFEO_DVECEL(qp_in->d+ii, j+nb[ii]+ng[ii]) -=
                        BLASFEO_DVECEL(qp_out->ux+ii, slack_index+nx[ii]+nu[ii]+ns[ii]);
            }
        }

        // NOTE: bounds on slacks can be skipped, since they are linear.
        // blasfeo_daxpy(2*ns[ii], -1.0, qp_out->ux+ii, nx[ii]+nu[ii], qp_in->d+ii, 2*nb[ii]+2*ng[ii], qp_in->d+ii, 2*nb[ii]+2*ng[ii]);

        // printf("SOC: qp_in->d final value\n");
        // blasfeo_print_exp_dvec(2*nb[ii]+2*ng[ii], qp_in->d+ii, 0);
    }

    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: SOC ocp_qp_in at iteration %d\n", sqp_iter);
        print_ocp_qp_in(qp_in);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
    ocp_nlp_sqp_dump_qp_in_to_file(qp_in, sqp_iter, 1);
#endif

    // solve QP
    // acados_tic(&timer1);
    int qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, qp_in, qp_out,
                                    opts->nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
    // NOTE: QP is not timed, since this computation time is attributed to globalization.

    // compute correct dual solution in case of Hessian regularization
    config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                        opts->nlp_opts->regularize, nlp_mem->regularize_mem);

    ocp_qp_out_get(qp_out, "qp_info", &qp_info_);
    int qp_iter = qp_info_->num_iter;

    // save statistics of last qp solver call
    // TODO: SOC QP solver call should be warm / hot started!
    if (sqp_iter+1 < mem->stat_m)
    {
        // mem->stat[mem->stat_n*(sqp_iter+1)+4] = qp_status;
        // add qp_iter; should maybe be in a seperate statistic
        mem->stat[mem->stat_n*(sqp_iter+1)+5] += qp_iter;
    }

    // compute external QP residuals (for debugging)
    if (opts->ext_qp_res)
    {
        ocp_qp_res_compute(qp_in, qp_out, work->qp_res, work->qp_res_ws);
        if (sqp_iter+1 < mem->stat_m)
            ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*(sqp_iter+1)+7));
    }

    if (nlp_opts->print_level > 3)
    {
        printf("\n\nSQP: SOC ocp_qp_out at iteration %d\n", sqp_iter);
        print_ocp_qp_out(qp_out);
    }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
        ocp_nlp_sqp_dump_qp_out_to_file(qp_out, sqp_iter, 1);
#endif

    // exit conditions on QP status
    if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
    {
#ifndef ACADOS_SILENT
        printf("\nQP solver returned error status %d in SQP iteration %d for SOC QP in QP iteration %d.\n",
            qp_status, sqp_iter, qp_iter);
#endif
        if (nlp_opts->print_level > 1)
        {
            printf("\nFailed to solve the following QP:\n");
            if (nlp_opts->print_level > 3)
                print_ocp_qp_in(qp_in);
        }

        mem->status = ACADOS_QP_FAILURE;
        mem->sqp_iter = sqp_iter;

        return ACADOS_QP_FAILURE;
    }
    return true;
}

static void ocp_nlp_sqp_reset_timers(ocp_nlp_sqp_memory *mem)
{
    mem->time_qp_sol = 0.0;
    mem->time_qp_solver_call = 0.0;
    mem->time_qp_xcond = 0.0;
    mem->time_lin = 0.0;
    mem->time_reg = 0.0;
    mem->time_glob = 0.0;
    mem->time_sim = 0.0;
    mem->time_sim_la = 0.0;
    mem->time_sim_ad = 0.0;
}

static double get_l1_infeasibility(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_sqp_memory *mem)
{
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    int i;
    int j;

    // compute current l1 infeasibility
    double tmp;
    struct blasfeo_dvec *tmp_fun_vec;
    double dyn_l1_infeasibility = 0.0;
    for(i=0; i<N; i++)
    {
        tmp_fun_vec = config->dynamics[i]->memory_get_fun_ptr(nlp_mem->dynamics[i]);
        for(j=0; j<nx[i+1]; j++)
        {
            dyn_l1_infeasibility += fabs(BLASFEO_DVECEL(tmp_fun_vec, j));
        }
    }

    double constraint_l1_infeasibility = 0.0;
    for(i=0; i<=N; i++)
    {
        tmp_fun_vec = config->constraints[i]->memory_get_fun_ptr(nlp_mem->constraints[i]);
        // tmp_fun_vec = out->t+i;
        for (j=0; j<2*ni[i]; j++)
        {
            tmp = BLASFEO_DVECEL(tmp_fun_vec, j);
            if (tmp > 0.0)
            {
                constraint_l1_infeasibility += tmp;
            }
        }
    }
    return dyn_l1_infeasibility + constraint_l1_infeasibility;
}

/************************************************
 * output functions
 ************************************************/
static void print_iteration_header(ocp_nlp_opts* opts){
    if(opts->globalization == FUNNEL_L1PEN_LINESEARCH)
    {
        printf("%6s | %11s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %12s | %10s | %10s | %10s | %10s\n",
        "iter.",
        "objective",
        "res_eq",
        "res_ineq",
        "res_stat",
        "res_comp",
        "alpha",
        "step_norm",
        "LM_reg.",
        "funnel width",
        "penalty",
        "qp_status",
        "qp_iter",
        "iter. type");
    }
    else
    {
        printf("# it\tstat\t\teq\t\tineq\t\tcomp\t\tqp_stat\tqp_iter\talpha\n");
    }
}

static void print_iteration(ocp_nlp_opts* opts,
                    double obj,
                    int iter_count,
                    double infeas_eq,
                    double infeas_ineq,
                    double stationarity,
                    double complementarity,
                    double alpha,
                    double step_norm,
                    double reg_param,
                    double funnel_width,
                    double penalty_parameter,
                    int qp_status,
                    int qp_iter,
                    char iter_type)
{
    if ((iter_count % 10 == 0)){
        print_iteration_header(opts);
    }
    if (opts->globalization == FUNNEL_L1PEN_LINESEARCH)
    {
        printf("%6i | %11.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %12.4e | %10.4e | %10i | %10i | %10c\n",
        iter_count,
        obj,
        infeas_eq,
        infeas_ineq,
        stationarity,
        complementarity,
        alpha,
        step_norm,
        reg_param,
        funnel_width,
        penalty_parameter,
        qp_status,
        qp_iter,
        iter_type);
    }
    else
    {
        printf("%i\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n",
        iter_count,
        stationarity,
        infeas_eq,
        infeas_ineq,
        complementarity,
        qp_status,
        qp_iter,
        alpha);
    }
}

/************************************************
 * termination criterion
 ************************************************/
static bool check_termination(int n_iter, ocp_nlp_res *nlp_res, ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_opts *opts)
{
    // check for nans
    if (isnan(nlp_res->inf_norm_res_stat) || isnan(nlp_res->inf_norm_res_eq) ||
            isnan(nlp_res->inf_norm_res_ineq) || isnan(nlp_res->inf_norm_res_comp))
    {
        mem->status = ACADOS_NAN_DETECTED;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Stopped: NaN detected in iterate.\n");
        }
        return true;
    }

    // check for maximum iterations
    if (!opts->eval_residual_at_max_iter && n_iter >= opts->max_iter)
    {
        mem->status = ACADOS_MAXITER;
        if (opts->nlp_opts->print_level > 0){
            printf("Stopped: Maximum Iterations Reached.\n");
        }
        return true;
    }

    // check if solved to tolerance
    if ((nlp_res->inf_norm_res_stat < opts->tol_stat) &&
        (nlp_res->inf_norm_res_eq < opts->tol_eq) &&
        (nlp_res->inf_norm_res_ineq < opts->tol_ineq) &&
        (nlp_res->inf_norm_res_comp < opts->tol_comp))
    {
        mem->status = ACADOS_SUCCESS;
        if (opts->nlp_opts->print_level > 0)
        {
            printf("Optimal Solution found! Converged to KKT point.\n");
        }
        return true;
    }

    // check for small step
    if (opts->tol_min_step_norm > 0.0 && (n_iter > 0) && (mem->step_norm < opts->tol_min_step_norm))
    {
        if (opts->nlp_opts->print_level > 0)
        {
            if (nlp_res->inf_norm_res_eq < opts->tol_eq && nlp_res->inf_norm_res_ineq < opts->tol_ineq)
            {
                printf("Stopped: Converged To Feasible Point. Step size is < tol_eq.\n");
            }
            else
            {
                printf("Stopped: Converged To Infeasible Point. Step size is < tol_eq.\n");
            }
        }
        mem->status = ACADOS_MINSTEP;
        return true;
    }

    // check for unbounded problem
    if (mem->nlp_mem->cost_value <= opts->tol_unbounded)
    {
        mem->status = ACADOS_UNBOUNDED;
        if (opts->nlp_opts->print_level > 0){
            printf("Stopped: Problem seems to be unbounded.\n");
        }
        return true;
    }

    // check for maximum iterations
    if (n_iter >= opts->max_iter)
    {
        mem->status = ACADOS_MAXITER;
        if (opts->nlp_opts->print_level > 0){
            printf("Stopped: Maximum Iterations Reached.\n");
        }
        return true;
    }

    return false;
}

/************************************************
 * funnel functions
 ************************************************/
static void debug_output(ocp_nlp_opts *opts, char* message, int print_level)
{
    if (opts->print_level > print_level)
    {
        printf("%s", message); //debugging output
    }
}
static void debug_output_double(ocp_nlp_opts *opts, char* message, double value, int print_level)
{
    if (opts->print_level > print_level)
    {
        printf("%s: %f\n", message, value); //debugging output
    }
}

static void initialize_funnel_width(ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_opts *opts, double initial_infeasibility)
{
    mem->funnel_width = fmax(opts->funnel_initialization_upper_bound,
                            opts->funnel_initialization_increase_factor*initial_infeasibility);
}

static void initialize_funnel_penalty_parameter(ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_opts *opts)
{
    mem->funnel_penalty_parameter = opts->funnel_initial_penalty_parameter;
}

static void update_funnel_penalty_parameter(ocp_nlp_sqp_memory *mem,
                                            ocp_nlp_sqp_opts *opts,
                                            double pred_f, double pred_h)
{
    if (mem->funnel_penalty_parameter * pred_f + pred_h < opts->funnel_penalty_eta * pred_h)
    {
        mem->funnel_penalty_parameter = fmin(opts->funnel_penalty_contraction * mem->funnel_penalty_parameter,
                                             ((1-opts->funnel_penalty_eta) * pred_h) / (-pred_f));
    }
    // else: do not decrease penalty parameter
}

static void decrease_funnel(ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_opts *opts, double trial_infeasibility, double current_infeasibility)
{
    mem->funnel_width = (1-opts->funnel_kappa) * trial_infeasibility + opts->funnel_kappa * mem->funnel_width;
}

static bool is_iterate_inside_of_funnel(ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_opts *opts, double infeasibility)
{
    if (infeasibility <= mem->funnel_width)
    {
        return true;
    }
    else
    {
        return false;
    }
}

static bool is_funnel_sufficient_decrease_satisfied(ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_opts *opts, double infeasibility)
{
    if (infeasibility <= opts->funnel_sufficient_decrease_factor* mem->funnel_width)
    {
        return true;
    }
    else
    {
        return false;
    }
}

static bool is_switching_condition_satisfied(ocp_nlp_sqp_opts *opts, double pred_optimality, double step_size, double pred_infeasibility)
{
    // if (step_size * pred_optimality >= opts->funnel_fraction_switching_condition * pred_infeasibility * pred_infeasibility)
    if (step_size * pred_optimality >= opts->funnel_fraction_switching_condition * pred_infeasibility)
    {
        return true;
    }
    else
    {
        return false;
    }
}

static bool is_f_type_armijo_condition_satisfied(ocp_nlp_sqp_opts *opts,
                                                    double negative_ared,
                                                    double pred,
                                                    double alpha)
{
    if (negative_ared <= fmin(-opts->nlp_opts->eps_sufficient_descent * alpha * fmax(pred, 0) + 1e-18, 0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

static bool is_trial_iterate_acceptable_to_funnel(ocp_nlp_sqp_memory *mem,
                                                  ocp_nlp_sqp_opts *opts,
                                                  double pred, double ared, double alpha,
                                                  double current_infeasibility,
                                                  double trial_infeasibility,
                                                  double current_objective,
                                                  double trial_objective,
                                                  double current_merit,
                                                  double trial_merit,
                                                  double pred_merit)
{
    bool accept_step = false;
    debug_output_double(opts->nlp_opts, "current objective", current_objective, 2);
    debug_output_double(opts->nlp_opts, "current infeasibility", current_infeasibility, 2); 
    debug_output_double(opts->nlp_opts, "trial objective", trial_objective, 2); 
    debug_output_double(opts->nlp_opts, "trial infeasibility", trial_infeasibility, 2); 
    debug_output_double(opts->nlp_opts, "pred", pred, 2); 

    if(is_iterate_inside_of_funnel(mem, opts, trial_infeasibility))
    {
        debug_output(opts->nlp_opts, "Trial iterate is INSIDE of funnel\n", 1);
        if (!mem->funnel_penalty_mode)
        {
            debug_output(opts->nlp_opts, "Penalty Mode not active!\n", 1);
            if (is_switching_condition_satisfied(opts, pred, alpha, current_infeasibility))
            {
                debug_output(opts->nlp_opts, "Switching condition IS satisfied!\n", 1);
                if (is_f_type_armijo_condition_satisfied(opts, -ared, pred, alpha))
                {
                    debug_output(opts->nlp_opts, "f-type step: Armijo condition satisfied\n", 1);
                    accept_step = true;
                    mem->funnel_iter_type = 'f';
                }
                else
                {
                    debug_output(opts->nlp_opts, "f-type step: Armijo condition NOT satisfied\n", 1);
                }

            }
            else if (is_funnel_sufficient_decrease_satisfied(mem, opts, trial_infeasibility))
            {
                debug_output(opts->nlp_opts, "Switching condition is NOT satisfied!\n", 1);
                debug_output(opts->nlp_opts, "h-type step: funnel suff. decrease satisfied!\n", 1);
                accept_step = true;
                mem->funnel_iter_type = 'h';
                decrease_funnel(mem, opts, trial_infeasibility, current_infeasibility);
            }
            else
            {
                debug_output(opts->nlp_opts, "Switching condition is NOT satisfied!\n", 1);
                debug_output(opts->nlp_opts, "Entered penalty check!\n", 1);
                //TODO move to function and test more
                if (trial_merit <= current_merit + opts->nlp_opts->eps_sufficient_descent * alpha * pred_merit)
                {
                    debug_output(opts->nlp_opts, "Penalty Function accepted\n", 1);
                    accept_step = true;
                    mem->funnel_iter_type = 'b';
                    mem->funnel_penalty_mode = true;
                }
            }
        }
        else
        {
            debug_output(opts->nlp_opts, "Penalty mode active\n", 1);
            if (trial_merit <= current_merit + opts->nlp_opts->eps_sufficient_descent * alpha * pred_merit)
            {
                debug_output(opts->nlp_opts, "p-type step: accepted iterate\n", 1);
                accept_step = true;
                mem->funnel_iter_type = 'p';

                if (is_funnel_sufficient_decrease_satisfied(mem, opts, trial_infeasibility))
                {
                    decrease_funnel(mem, opts, trial_infeasibility, current_infeasibility);
                    mem->funnel_penalty_mode = false;
                }
            }
        }
    }
    else
    {
        debug_output(opts->nlp_opts, "Trial iterate is NOT INSIDE of funnel\n", 1);
    }
    return accept_step;
}

static int ocp_nlp_sqp_backtracking_line_search(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *mem_, void *work_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;

    int N = dims->N;
    double pred = -nlp_mem->qp_cost_value;
    double pred_merit = 0.0; // Calculate this here
    double alpha = 1.0;
    double trial_cost;
    double trial_infeasibility = 0.0;
    double ared;
    bool accept_step;
    double current_infeasibility = mem->l1_infeasibility;
    double current_cost = nlp_mem->cost_value;
    double current_merit = mem->funnel_penalty_parameter*current_cost + current_infeasibility;

    // do the penalty parameter update here .... might be changed later
    update_funnel_penalty_parameter(mem, opts, pred, mem->l1_infeasibility);

    int i;

    while (true)
    {
        // Calculate trial iterate: trial_iterate = current_iterate + alpha * direction
        ocp_nlp_update_variables_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem,
                                     nlp_work, nlp_work->tmp_nlp_out, alpha);

        ///////////////////////////////////////////////////////////////////////
        // Evaluate cost function at trial iterate
        // set evaluation point to tmp_nlp_out
        ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, nlp_in, nlp_work->tmp_nlp_out, nlp_mem);
        // compute trial dynamics value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
        for (i=0; i<N; i++)
        {
            // dynamics: Note has to be first, because cost_integration might be used.
            config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i], nlp_in->dynamics[i],
                                            nlp_opts->dynamics[i], nlp_mem->dynamics[i], nlp_work->dynamics[i]);
        }
        // compute trial objective function value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
        for (i=0; i<=N; i++)
        {
            // cost
            config->cost[i]->compute_fun(config->cost[i], dims->cost[i], nlp_in->cost[i], nlp_opts->cost[i],
                                        nlp_mem->cost[i], nlp_work->cost[i]);
        }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
        for (i=0; i<=N; i++)
        {
            // constr
            config->constraints[i]->compute_fun(config->constraints[i], dims->constraints[i],
                                                nlp_in->constraints[i], nlp_opts->constraints[i],
                                                nlp_mem->constraints[i], nlp_work->constraints[i]);
        }
        // reset evaluation point to SQP iterate
        ocp_nlp_set_primal_variable_pointers_in_submodules(config, dims, nlp_in, nlp_out, nlp_mem);

        double *tmp_fun;
        // Calculate the trial objective and constraint violation
        trial_cost = 0.0;
        for(i=0; i<=N; i++)
        {
            tmp_fun = config->cost[i]->memory_get_fun_ptr(nlp_mem->cost[i]);
            trial_cost += *tmp_fun;
        }
        trial_infeasibility = get_l1_infeasibility(config, dims, mem);

        ///////////////////////////////////////////////////////////////////////
        // Evaluate merit function at trial point
        double trial_merit = mem->funnel_penalty_parameter*trial_cost + trial_infeasibility;
        pred_merit = mem->funnel_penalty_parameter * pred + current_infeasibility;
        ared = nlp_mem->cost_value - trial_cost;

        // Funnel globalization
        accept_step = is_trial_iterate_acceptable_to_funnel(mem, opts,
                                                            pred, ared,
                                                            alpha, current_infeasibility,
                                                            trial_infeasibility, current_cost,
                                                            trial_cost, current_merit, trial_merit,
                                                            pred_merit);

        if (accept_step)
        {
            mem->alpha = alpha;
            nlp_mem->cost_value = trial_cost;
            mem->l1_infeasibility = trial_infeasibility;
            return 1;
        }

        if (alpha < opts->nlp_opts->alpha_min)
        {
            printf("Linesearch: Step size gets too small. Should enter penalty phase. \n");
            exit(1);
        }

        alpha *= opts->nlp_opts->alpha_reduction;

    }
}

/************************************************
 * functions
 ************************************************/

// MAIN OPTIMIZATION ROUTINE
int ocp_nlp_sqp(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    acados_timer timer0, timer1;
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_res *nlp_res = nlp_mem->nlp_res;

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_qp_in *qp_in = nlp_mem->qp_in;
    ocp_qp_out *qp_out = nlp_mem->qp_out;

    // zero timers
    double tmp_time;
    ocp_nlp_sqp_reset_timers(mem);

    int N = dims->N;
    int ii;
    int qp_status = 0;
    int qp_iter = 0;
    mem->alpha = 0.0;
    mem->funnel_iter_type = '-';
    mem->status = ACADOS_SUCCESS;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
#endif

    ocp_nlp_initialize_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    /************************************************
     * main sqp loop
     ************************************************/
    int sqp_iter = 0;
    double reg_param_memory = 0.0;
    double funnel_width_memory = 0.0;
    double funnel_penalty_param_memory = opts->funnel_initial_penalty_parameter;
    initialize_funnel_penalty_parameter(mem, opts);
    if (nlp_opts->globalization == FUNNEL_L1PEN_LINESEARCH)
    {
        printf("Note: The funnel globalization is still under development.\n");
        printf("If you encouter problems or bugs, please report to the acados developers!\n");
    }

    for (; sqp_iter <= opts->max_iter; sqp_iter++) // <= needed such that after last iteration KKT residuals are checked before max_iter is thrown.
    {
        // We always evaluate the residuals until the last iteration
        // If the option "eval_residual_at_max_iter" is set, we also
        // evaluate the residuals after the last iteration.
        if (sqp_iter != opts->max_iter || opts->eval_residual_at_max_iter)
        {
            /* Prepare the QP data */
            // linearize NLP and update QP matrices
            acados_tic(&timer1);
            ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            if (nlp_opts->with_adaptive_levenberg_marquardt || nlp_opts->globalization != FIXED_STEP)
            {
                ocp_nlp_get_cost_value_from_submodules(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
            }
            ocp_nlp_add_levenberg_marquardt_term(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, mem->alpha, sqp_iter);

            mem->time_lin += acados_toc(&timer1);

            // get timings from integrator
            for (ii=0; ii<N; ii++)
            {
                config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], "time_sim", &tmp_time);
                mem->time_sim += tmp_time;
                config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], "time_sim_la", &tmp_time);
                mem->time_sim_la += tmp_time;
                config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], "time_sim_ad", &tmp_time);
                mem->time_sim_ad += tmp_time;
            }

            // update QP rhs for SQP (step prim var, abs dual var)
            ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

            // compute nlp residuals
            ocp_nlp_res_compute(dims, nlp_in, nlp_out, nlp_res, nlp_mem);
            ocp_nlp_res_get_inf_norm(nlp_res, &nlp_out->inf_norm_res);

            if (nlp_opts->globalization == FUNNEL_L1PEN_LINESEARCH && sqp_iter == 0)
            {
                mem->l1_infeasibility = get_l1_infeasibility(config, dims, mem);
            }
        }

        // initialize funnel if FUNNEL_L1PEN_LINESEARCH used
        if (sqp_iter == 0 && nlp_opts->globalization == FUNNEL_L1PEN_LINESEARCH){
            initialize_funnel_width(mem, opts, mem->l1_infeasibility);
        }
        funnel_width_memory = mem->funnel_width;
        funnel_penalty_param_memory = mem->funnel_penalty_parameter;

        // save statistics
        if (sqp_iter < mem->stat_m)
        {
            mem->stat[mem->stat_n*sqp_iter+0] = nlp_res->inf_norm_res_stat;
            mem->stat[mem->stat_n*sqp_iter+1] = nlp_res->inf_norm_res_eq;
            mem->stat[mem->stat_n*sqp_iter+2] = nlp_res->inf_norm_res_ineq;
            mem->stat[mem->stat_n*sqp_iter+3] = nlp_res->inf_norm_res_comp;
        }

        // Output
        if (nlp_opts->print_level > 0)
        {
            print_iteration(nlp_opts, nlp_mem->cost_value, sqp_iter, nlp_res->inf_norm_res_eq,
                            nlp_res->inf_norm_res_ineq, nlp_res->inf_norm_res_stat,
                            nlp_res->inf_norm_res_comp, mem->alpha, mem->step_norm,
                            reg_param_memory, funnel_width_memory, funnel_penalty_param_memory, qp_status,
                            qp_iter, mem->funnel_iter_type);
        }
        reg_param_memory = nlp_opts->levenberg_marquardt;

        // Termination
        if (check_termination(sqp_iter, nlp_res, mem, opts))
        {
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            mem->sqp_iter = sqp_iter;
            mem->time_tot = acados_toc(&timer0);
            return mem->status;
        }


        // regularize Hessian
        acados_tic(&timer1);
        config->regularize->regularize(config->regularize, dims->regularize,
                                               nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        /* solve QP */
        // (typically) no warm start at first iteration
        if (sqp_iter == 0 && !opts->warm_start_first_qp)
        {
            int tmp_int = 0;
            config->qp_solver->opts_set(config->qp_solver, nlp_opts->qp_solver_opts,
                                         "warm_start", &tmp_int);
        }
        // Show input to QP
        if (nlp_opts->print_level > 3)
        {
            printf("\n\nSQP: ocp_qp_in at iteration %d\n", sqp_iter);
            print_ocp_qp_in(qp_in);
        }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
        ocp_nlp_sqp_dump_qp_in_to_file(qp_in, sqp_iter, 0);
#endif
        // solve qp
        acados_tic(&timer1);
        qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, qp_in, qp_out,
                                        nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
        mem->time_qp_sol += acados_toc(&timer1);

        qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_solver_call", &tmp_time);
        mem->time_qp_solver_call += tmp_time;
        qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_xcond", &tmp_time);
        mem->time_qp_xcond += tmp_time;

        // compute correct dual solution in case of Hessian regularization
        acados_tic(&timer1);
        config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                             nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        // restore default warm start
        if (sqp_iter==0)
        {
            config->qp_solver->opts_set(config->qp_solver, nlp_opts->qp_solver_opts,
                                        "warm_start", &opts->qp_warm_start);
        }

        if (nlp_opts->print_level > 3)
        {
            printf("\n\nSQP: ocp_qp_out at iteration %d\n", sqp_iter);
            print_ocp_qp_out(qp_out);
        }

#if defined(ACADOS_DEBUG_SQP_PRINT_QPS_TO_FILE)
        ocp_nlp_sqp_dump_qp_out_to_file(qp_out, sqp_iter, 0);
#endif

        qp_info *qp_info_;
        ocp_qp_out_get(qp_out, "qp_info", &qp_info_);
        qp_iter = qp_info_->num_iter;

        // save statistics of last qp solver call
        if (sqp_iter+1 < mem->stat_m)
        {
            mem->stat[mem->stat_n*(sqp_iter+1)+4] = qp_status;
            mem->stat[mem->stat_n*(sqp_iter+1)+5] = qp_iter;
        }

        // compute external QP residuals (for debugging)
        if (opts->ext_qp_res)
        {
            ocp_qp_res_compute(qp_in, qp_out, work->qp_res, work->qp_res_ws);
            if (sqp_iter+1 < mem->stat_m)
                ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*(sqp_iter+1)+7));
        }

        // exit conditions on QP status
        if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
        {
            if (nlp_opts->print_level > 0)
            {
                printf("%i\t%e\t%e\t%e\t%e.\n", sqp_iter, nlp_res->inf_norm_res_stat,
                    nlp_res->inf_norm_res_eq, nlp_res->inf_norm_res_ineq,
                    nlp_res->inf_norm_res_comp );
                printf("\n\n");
            }
            // increment sqp_iter to return full statistics and improve output below.
            sqp_iter++;

#ifndef ACADOS_SILENT
            printf("\nQP solver returned error status %d in SQP iteration %d, QP iteration %d.\n",
                   qp_status, sqp_iter, qp_iter);
#endif
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            if (nlp_opts->print_level > 1)
            {
                printf("\n Failed to solve the following QP:\n");
                if (nlp_opts->print_level)
                    print_ocp_qp_in(qp_in);
            }

            mem->status = ACADOS_QP_FAILURE;
            mem->sqp_iter = sqp_iter;
            mem->time_tot = acados_toc(&timer0);

            return mem->status;
        }

        // Calculate optimal QP objective (needed for globalization)
        if (nlp_opts->globalization == FUNNEL_L1PEN_LINESEARCH)
        {
            nlp_mem->qp_cost_value = ocp_nlp_sqp_compute_qp_objective_value(dims, qp_in, qp_out, nlp_work, nlp_mem, opts);
        }

        // Compute the step norm
        if (opts->tol_min_step_norm > 0.0 || nlp_opts->log_primal_step_norm)
        {
            mem->step_norm = ocp_qp_out_compute_primal_nrm_inf(nlp_mem->qp_out);
            if (nlp_opts->log_primal_step_norm)
                mem->primal_step_norm[sqp_iter] = nlp_mem->qp_cost_value;
        }
        /* end solve QP */

        /* globalization */
        // NOTE on timings: currently all within globalization is accounted for within time_glob.
        //   QP solver times could be also attributed there alternatively. Cleanest would be to save them seperately.
        acados_tic(&timer1);
        if (nlp_opts->globalization == FUNNEL_L1PEN_LINESEARCH)
        {
            bool linesearch_success = 1;
            linesearch_success = ocp_nlp_sqp_backtracking_line_search(config, dims, nlp_in, nlp_out, mem, work, opts);
             // Copy new iterate to nlp_out
            if (linesearch_success)
            {
                // in case line search fails, we do not want to copy trial iterates!
                copy_ocp_nlp_out(dims, work->nlp_work->tmp_nlp_out, nlp_out);
            }
            mem->time_glob += acados_toc(&timer1);
        }
        else
        {
            bool do_line_search = true;
            if (nlp_opts->globalization_use_SOC && nlp_opts->globalization == MERIT_BACKTRACKING)
            {
                do_line_search = ocp_nlp_soc_line_search(config, dims, nlp_in, nlp_out, opts, mem, work, sqp_iter);
                if (mem->status == ACADOS_QP_FAILURE)
                {
#if defined(ACADOS_WITH_OPENMP)
                    // restore number of threads
                    omp_set_num_threads(num_threads_bkp);
#endif
                    mem->time_tot = acados_toc(&timer0);
                    return mem->status;
                }
            }

            if (do_line_search)
            {
                mem->alpha = ocp_nlp_line_search(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, 0, sqp_iter);
            }
            mem->time_glob += acados_toc(&timer1);
            mem->stat[mem->stat_n*(sqp_iter+1)+6] = mem->alpha;

            // update variables
            ocp_nlp_update_variables_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, nlp_out, mem->alpha);
        }

    }  // end SQP loop

    if (nlp_opts->print_level > 0)
    {
        printf("Warning: The solver should never reach this part of the function!\n");
    }
    return mem->status;
}

double ocp_nlp_sqp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                ocp_nlp_workspace *nlp_work, ocp_nlp_memory *nlp_mem, ocp_nlp_sqp_opts *opts)
{
    // Compute the QP objective function value
    double qp_cost = 0.0;
    int i, nux, ns;
    int N = dims->N;
    // Sum over stages 0 to N
    for (i = 0; i <= N; i++)
    {
        nux = dims->nx[i] + dims->nu[i];
        ns = dims->ns[i];
        if (opts->funnel_type_switching_condition)
        {
            // Calculate 0.5 * d.T H d
            blasfeo_dsymv_l(nux, 0.5, &qp_in->RSQrq[i], 0, 0, &qp_out->ux[i], 0,
                            0.0, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);
            qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);

            // slack QP objective value, compare to computation in cost modules;
            // tmp_nv = 2 * z + Z .* slack;
            blasfeo_dveccpsc(2*ns, 2.0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
            blasfeo_dvecmulacc(2*ns, &qp_in->Z[i], 0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
            // qp_cost += .5 * (tmp_nv .* slack)
            qp_cost += 0.5 * blasfeo_ddot(2*ns, &nlp_work->tmp_nv, 0, &qp_out->ux[i], nux);
        }
        // Calculate g.T d
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &qp_in->rqz[i], 0);

        // Calculate gradient of slacks
        qp_cost += blasfeo_ddot(2 * ns, &qp_out->ux[i], nux, &qp_in->rqz[i], nux);
    }
    return qp_cost;
}

void ocp_nlp_sqp_memory_reset_qp_solver(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
    void *opts_, void *mem_, void *work_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    // printf("in ocp_nlp_sqp_memory_reset_qp_solver\n\n");
    config->qp_solver->memory_reset(qp_solver, dims->qp_solver,
        nlp_mem->qp_in, nlp_mem->qp_out, opts->nlp_opts->qp_solver_opts,
        nlp_mem->qp_solver_mem, nlp_work->qp_work);
}


int ocp_nlp_sqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    nlp_mem->workspace_size = ocp_nlp_workspace_calculate_size(config, dims, opts->nlp_opts);

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_sqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    return ocp_nlp_precompute_common(config, dims, nlp_in, nlp_out, opts->nlp_opts, nlp_mem, nlp_work);
}




void ocp_nlp_sqp_eval_param_sens(void *config_, void *dims_, void *opts_, void *mem_, void *work_,
                                 char *field, int stage, int index, void *sens_nlp_out_)
{
    acados_timer timer0;
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_out *sens_nlp_out = sens_nlp_out_;

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_common_eval_param_sens(config, dims, opts->nlp_opts, nlp_mem, nlp_work,
                                 field, stage, index, sens_nlp_out);

    mem->time_solution_sensitivities = acados_toc(&timer0);

    return;
}


void ocp_nlp_sqp_eval_lagr_grad_p(void *config_, void *dims_, void *nlp_in_, void *opts_, void *mem_, void *work_,
                                 const char *field, void *grad_p)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_in *nlp_in = nlp_in_;

    ocp_nlp_sqp_workspace *work = work_;
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    ocp_nlp_common_eval_lagr_grad_p(config, dims, nlp_in, opts->nlp_opts, nlp_mem, nlp_work,
                                 field, grad_p);

    return;
}


void ocp_nlp_sqp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_sqp_memory *mem = mem_;


    if (!strcmp("sqp_iter", field) || !strcmp("nlp_iter", field))
    {
        int *value = return_value_;
        *value = mem->sqp_iter;
    }
    else if (!strcmp("status", field))
    {
        int *value = return_value_;
        *value = mem->status;
    }
    else if (!strcmp("time_tot", field) || !strcmp("tot_time", field))
    {
        double *value = return_value_;
        *value = mem->time_tot;
    }
    else if (!strcmp("time_qp_sol", field) || !strcmp("time_qp", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_sol;
    }
    else if (!strcmp("time_qp_solver", field) || !strcmp("time_qp_solver_call", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_solver_call;
    }
    else if (!strcmp("time_qp_xcond", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_xcond;
    }
    else if (!strcmp("time_lin", field))
    {
        double *value = return_value_;
        *value = mem->time_lin;
    }
    else if (!strcmp("time_reg", field))
    {
        double *value = return_value_;
        *value = mem->time_reg;
    }
    else if (!strcmp("time_glob", field))
    {
        double *value = return_value_;
        *value = mem->time_glob;
    }
    else if (!strcmp("time_solution_sensitivities", field))
    {
        double *value = return_value_;
        *value = mem->time_solution_sensitivities;
    }
    else if (!strcmp("time_sim", field))
    {
        double *value = return_value_;
        *value = mem->time_sim;
    }
    else if (!strcmp("time_sim_la", field))
    {
        double *value = return_value_;
        *value = mem->time_sim_la;
    }
    else if (!strcmp("time_sim_ad", field))
    {
        double *value = return_value_;
        *value = mem->time_sim_ad;
    }
    else if (!strcmp("time_preparation", field))
    {
        double *value = return_value_;
        *value = 0.0;
    }
    else if (!strcmp("time_feedback", field))
    {
        double *value = return_value_;
        *value = mem->time_tot;
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("primal_step_norm", field))
    {
        if (mem->primal_step_norm == NULL)
        {
            printf("\nerror: options log_primal_step_norm was not set\n");
            exit(1);
        }
        else
        {
            double *value = return_value_;
            for (int ii=0; ii<mem->sqp_iter; ii++)
            {
                value[ii] = mem->primal_step_norm[ii];
            }
        }
    }

    else if (!strcmp("statistics", field))
    {
        int n_row = mem->stat_m<mem->sqp_iter+1 ? mem->stat_m : mem->sqp_iter+1;
        double *value = return_value_;
        for (int ii=0; ii<n_row; ii++)
        {
            value[ii+0] = ii;
            for (int jj=0; jj<mem->stat_n; jj++)
                value[ii+(jj+1)*n_row] = mem->stat[jj+ii*mem->stat_n];
        }
    }
    else if (!strcmp("stat_m", field))
    {
        int *value = return_value_;
        *value = mem->stat_m;
    }
    else if (!strcmp("stat_n", field))
    {
        int *value = return_value_;
        *value = mem->stat_n;
    }
    else if (!strcmp("nlp_mem", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem;
    }
    else if (!strcmp("qp_xcond_dims", field))
    {
        void **value = return_value_;
        *value = dims->qp_solver->xcond_dims;
    }
    else if (!strcmp("nlp_res", field))
    {
        ocp_nlp_res **value = return_value_;
        *value = mem->nlp_mem->nlp_res;
    }
    else if (!strcmp("qp_xcond_in", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_solver_mem->xcond_qp_in;
    }
    else if (!strcmp("qp_xcond_out", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_solver_mem->xcond_qp_out;
    }
    else if (!strcmp("qp_in", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_in;
    }
    else if (!strcmp("qp_out", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_out;
    }
    else if (!strcmp("qp_iter", field))
    {
        config->qp_solver->memory_get(config->qp_solver,
            mem->nlp_mem->qp_solver_mem, "iter", return_value_);
    }
    else if (!strcmp("qp_status", field))
    {
        config->qp_solver->memory_get(config->qp_solver,
            mem->nlp_mem->qp_solver_mem, "status", return_value_);
    }
    else if (!strcmp("res_stat", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_stat;
    }
    else if (!strcmp("res_eq", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_eq;
    }
    else if (!strcmp("res_ineq", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_ineq;
    }
    else if (!strcmp("res_comp", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_comp;
    }
    else if (!strcmp("cost_value", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->cost_value;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_get\n", field);
        exit(1);
    }
}



void ocp_nlp_sqp_opts_get(void *config_, void *dims_, void *opts_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;

    if (!strcmp("nlp_opts", field))
    {
        void **value = return_value_;
        *value = opts->nlp_opts;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_opts_get\n", field);
        exit(1);
    }
}


void ocp_nlp_sqp_work_get(void *config_, void *dims_, void *work_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_sqp_workspace *work = work_;

    if (!strcmp("nlp_work", field))
    {
        void **value = return_value_;
        *value = work->nlp_work;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_sqp_work_get\n", field);
        exit(1);
    }
}



void ocp_nlp_sqp_terminate(void *config_, void *mem_, void *work_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_sqp_workspace *work = work_;

    config->qp_solver->terminate(config->qp_solver, mem->nlp_mem->qp_solver_mem, work->nlp_work->qp_work);
}


void ocp_nlp_sqp_config_initialize_default(void *config_)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    config->opts_calculate_size = &ocp_nlp_sqp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_sqp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_sqp_opts_initialize_default;
    config->opts_update = &ocp_nlp_sqp_opts_update;
    config->opts_set = &ocp_nlp_sqp_opts_set;
    config->opts_set_at_stage = &ocp_nlp_sqp_opts_set_at_stage;
    config->memory_calculate_size = &ocp_nlp_sqp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_sqp_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_sqp_workspace_calculate_size;
    config->evaluate = &ocp_nlp_sqp;
    config->memory_reset_qp_solver = &ocp_nlp_sqp_memory_reset_qp_solver;
    config->eval_param_sens = &ocp_nlp_sqp_eval_param_sens;
    config->eval_lagr_grad_p = &ocp_nlp_sqp_eval_lagr_grad_p;
    config->config_initialize_default = &ocp_nlp_sqp_config_initialize_default;
    config->precompute = &ocp_nlp_sqp_precompute;
    config->get = &ocp_nlp_sqp_get;
    config->opts_get = &ocp_nlp_sqp_opts_get;
    config->work_get = &ocp_nlp_sqp_work_get;
    config->terminate = &ocp_nlp_sqp_terminate;

    return;
}


// ??? @rien
//        for (int_t i = 0; i < N; i++)
//        {
//   ocp_nlp_dynamics_opts *dynamics_opts = opts->dynamics[i];
//            sim_opts *opts = dynamics_opts->sim_solver;
//            if (opts->scheme == NULL)
//                continue;
//            opts->sens_adj = (opts->scheme->type != exact);
//            if (nlp_in->freezeSens) {
//                // freeze inexact sensitivities after first SQP iteration !!
//                opts->scheme->freeze = true;
//            }
//        }
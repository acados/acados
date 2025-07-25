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
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"

#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_fixed_step.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"

/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_globalization_fixed_step_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_fixed_step_opts);

    size += ocp_nlp_globalization_opts_calculate_size(config, dims);

    return size;
}

void ocp_nlp_globalization_fixed_step_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_globalization_fixed_step_opts *opts = opts_;
    ocp_nlp_globalization_opts *globalization_opts = opts->globalization_opts;
    ocp_nlp_globalization_config *config = config_;

    opts->step_length = 1.0;
    ocp_nlp_globalization_opts_initialize_default(config, dims, globalization_opts);
    return;
}


void ocp_nlp_globalization_fixed_step_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_globalization_fixed_step_opts *opts = opts_;
    ocp_nlp_globalization_config *config = config_;

    if (!strcmp(field, "fixed_step_length"))
    {
        double* step_step = (double *) value;
        if (*step_step < 0.0 || *step_step > 1.0)
        {
            printf("\nerror: ocp_nlp_globalization_fixed_step_opts_set: invalid value for field fixed_step_length, need double in [0,1], got %f.", *step_step);
            exit(1);
        }
        opts->step_length = *step_step;
    }
    else
    {
        ocp_nlp_globalization_opts_set(config, opts->globalization_opts, field, value);
    }
    return;
}

void *ocp_nlp_globalization_fixed_step_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_globalization_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_globalization_fixed_step_opts *opts = (ocp_nlp_globalization_fixed_step_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_fixed_step_opts);

    opts->globalization_opts = ocp_nlp_globalization_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_globalization_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_globalization_fixed_step_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_globalization_fixed_step_memory_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_fixed_step_memory);

    return size;
}

void *ocp_nlp_globalization_fixed_step_memory_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_globalization_fixed_step_memory *mem = (ocp_nlp_globalization_fixed_step_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_fixed_step_memory);

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_globalization_fixed_step_memory_calculate_size(config_, dims_) >= c_ptr);

    return mem;
}


/************************************************
 * fixed step functions
 ************************************************/
int ocp_nlp_globalization_fixed_step_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *solver_mem, void *nlp_work_, void *nlp_opts_, double *step_size)
{
    ocp_nlp_config *config = nlp_config_;
    ocp_nlp_dims *dims = nlp_dims_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = nlp_mem_;
    ocp_nlp_workspace *nlp_work = nlp_work_;
    ocp_nlp_opts *nlp_opts = nlp_opts_;
    ocp_nlp_globalization_fixed_step_opts *opts = nlp_opts->globalization;

    ocp_qp_out *qp_out = nlp_mem->qp_out;
    double alpha = opts->step_length;

    if (nlp_opts->with_anderson_acceleration)
    {
        // convert qp_out to delta primal-dual step
        ocp_nlp_convert_primaldelta_absdual_step_to_delta_step(config, dims, nlp_out, qp_out);
        if (nlp_mem->iter == 0)
        {
            // store in anderson_step, prev_qp_out
            ocp_qp_out_copy(qp_out, nlp_mem->anderson_step);
            // update variables (TODO: DDP primals are different)
            ocp_nlp_update_variables_sqp_delta_primal_dual(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, alpha, nlp_mem->anderson_step);
        }
        else
        {
            // tmp_qp_out = d_{k+1} - d_k: qp_step - prev_qp_out
            ocp_qp_out_axpy(-1.0, nlp_mem->prev_qp_out, qp_out, nlp_work->tmp_qp_out);
            // compute gamma
            double gamma = ocp_nlp_compute_anderson_gamma(nlp_work, qp_out, nlp_work->tmp_qp_out);
            /* update anderson_step */
            // anderson_step *= -gamma
            ocp_qp_out_sc(-gamma, nlp_mem->anderson_step);
            // anderson_step += alpha * gamma * prev_qp_out
            ocp_qp_out_add(gamma*alpha, nlp_mem->prev_qp_out, nlp_mem->anderson_step);
            // anderson_step += (alpha - alpha * gamma) * qp_out
            ocp_qp_out_add(alpha-gamma*alpha, qp_out, nlp_mem->anderson_step);
            // update variables (TODO: DDP primals are different)
            ocp_nlp_update_variables_sqp_delta_primal_dual(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, alpha, nlp_mem->anderson_step);
        }
        // store prev qp step
        ocp_qp_out_copy(qp_out, nlp_mem->prev_qp_out);
        // step norm
        if (nlp_opts->log_primal_step_norm)
        {
            nlp_mem->primal_step_norm[nlp_mem->iter] = ocp_qp_out_compute_primal_nrm_inf(nlp_mem->anderson_step);
        }
        if (nlp_opts->log_dual_step_norm)
        {
            nlp_mem->dual_step_norm[nlp_mem->iter] = ocp_qp_out_compute_dual_nrm_inf(nlp_mem->anderson_step);
        }
    }
    else
    {
        config->step_update(config, dims, nlp_in, nlp_out, qp_out, nlp_opts, nlp_mem, nlp_work, nlp_out, solver_mem, alpha, opts->globalization_opts->full_step_dual);
    }
    *step_size = alpha;

    return ACADOS_SUCCESS;
}

void ocp_nlp_globalization_fixed_step_print_iteration_header()
{
    printf("%8s   ", "alpha");
}

void ocp_nlp_globalization_fixed_step_print_iteration(double objective_value,
                                                void* nlp_opts_,
                                                void* mem_)
{
    ocp_nlp_globalization_fixed_step_opts *opts = nlp_opts_;
    printf("%8.2e    ", opts->step_length);
}

int ocp_nlp_globalization_fixed_step_needs_objective_value()
{
    return 0;
}

int ocp_nlp_globalization_fixed_step_needs_qp_objective_value()
{
    return 0;
}

void ocp_nlp_globalization_fixed_step_initialize_memory(void *config_, void *dims_, void *nlp_mem_, void *nlp_opts_)
{
    return;
}

void ocp_nlp_globalization_fixed_step_config_initialize_default(ocp_nlp_globalization_config *config)
{
    // opts
    config->opts_calculate_size = &ocp_nlp_globalization_fixed_step_opts_calculate_size;
    config->opts_assign = &ocp_nlp_globalization_fixed_step_opts_assign;
    config->opts_initialize_default = &ocp_nlp_globalization_fixed_step_opts_initialize_default;
    config->opts_set = &ocp_nlp_globalization_fixed_step_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_globalization_fixed_step_memory_calculate_size;
    config->memory_assign = &ocp_nlp_globalization_fixed_step_memory_assign;
    // functions
    config->find_acceptable_iterate = &ocp_nlp_globalization_fixed_step_find_acceptable_iterate;
    config->print_iteration_header = &ocp_nlp_globalization_fixed_step_print_iteration_header;
    config->print_iteration = &ocp_nlp_globalization_fixed_step_print_iteration;
    config->needs_objective_value = &ocp_nlp_globalization_fixed_step_needs_objective_value;
    config->needs_qp_objective_value = &ocp_nlp_globalization_fixed_step_needs_qp_objective_value;
    config->initialize_memory = &ocp_nlp_globalization_fixed_step_initialize_memory;
}

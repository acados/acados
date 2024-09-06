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


#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_fixed_step.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"

// TODO: copy boilerblate..
// fix imports

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
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"


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
    ocp_nlp_globalization_opts_initialize_default(config_, dims_, opts_);

    return;
}


void ocp_nlp_globalization_fixed_step_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_globalization_fixed_step_opts *opts = opts_;
    ocp_nlp_globalization_config *config = config_;

    config->opts_set(config, opts->globalization_opts, field, value);

    return;
}

void *ocp_nlp_globalization_fixed_step_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_globalization_fixed_step_opts *opts = (ocp_nlp_globalization_fixed_step_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_fixed_step_opts);

    assert((char *) raw_memory + ocp_nlp_globalization_fixed_step_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_globalization_fixed_step_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_globalization_fixed_step_memory);

    return size;
}

void *ocp_nlp_globalization_fixed_step_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_globalization_fixed_step_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_globalization_fixed_step_memory *mem = (ocp_nlp_globalization_fixed_step_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_fixed_step_memory);

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_globalization_fixed_step_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}



void ocp_nlp_globalization_fixed_step_print_iteration_header()
{
    printf("# it\tstat\t\teq\t\tineq\t\tcomp\t\tqp_stat\tqp_iter\talpha\n");
}

// TODO: unified signature:
// -> move everything around.
// 1. residual_iter
// 2. int iter count
// 3. alpha etc. move to glob_memory. (void *)
void ocp_nlp_globalization_fixed_step_print_iteration(ocp_nlp_opts* opts,
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
        ocp_nlp_globalization_fixed_step_print_iteration_header(opts);
    }
    printf("%i\t%e\t%e\t%e\t%e\t%d\t%d\t%e\n",
        iter_count,
        stationarity,
        infeas_eq,        infeas_ineq,
        complementarity,
        qp_status,
        qp_iter,
        alpha);
}

int ocp_nlp_globalization_fixed_step_needs_objective_value()
{
    return 0;
}

int ocp_nlp_globalization_fixed_step_find_acceptable_iterate(void *nlp_config_, void *nlp_dims_, void *nlp_in_, void *nlp_out_, void *nlp_mem_, void *nlp_work_, void *nlp_opts_)
{
    ocp_nlp_config *nlp_config = nlp_config_;
    ocp_nlp_dims *nlp_dims = nlp_dims_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = nlp_mem_;
    ocp_nlp_workspace *nlp_work = nlp_work_;
    ocp_nlp_opts *nlp_opts = nlp_opts_;
    
    int sqp_iter = 1;
    bool do_line_search = true;
//     if (nlp_opts->globalization->globalization_use_SOC)
//     {
//         do_line_search = ocp_nlp_soc_line_search(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, sqp_iter);
// //         if (nlp_mem->status == ACADOS_QP_FAILURE)
// //         {
// // #if defined(ACADOS_WITH_OPENMP)
// //             // restore number of threads
// //             omp_set_num_threads(num_threads_bkp);
// // #endif
// //             // mem->time_tot = acados_toc(&timer0);
// //             return mem->status;
// //         }
//     }

    if (do_line_search)
    {
        int line_search_status;
        line_search_status = ocp_nlp_line_search(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, sqp_iter, &mem->alpha);
        if (line_search_status == ACADOS_NAN_DETECTED)
        {
            mem->status = ACADOS_NAN_DETECTED;
            return mem->status;
        }
    }
    // mem->time_glob += acados_toc(&timer1);
    // nlp_mem->stat[mem->stat_n*(sqp_iter+1)+6] = mem->alpha;

    // update variables
    ocp_nlp_update_variables_sqp(nlp_config, nlp_dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, nlp_out, mem->alpha);
}

void ocp_nlp_globalization_fixed_step_config_initialize_default(ocp_nlp_globalization_config *config)
{
    // opts
    config->opts_calculate_size = &ocp_nlp_globalization_fixed_step_opts_calculate_size;
    config->opts_assign = &ocp_nlp_globalization_fixed_step_opts_assign;
    config->opts_initialize_default = &ocp_nlp_globalization_fixed_step_opts_initialize_default;
    config->opts_set = &ocp_nlp_globalization_fixed_step_opts_set;
    // functions
    config->find_acceptable_iterate = &ocp_nlp_globalization_fixed_step_find_acceptable_iterate;
    config->print_iteration_header = &ocp_nlp_globalization_fixed_step_print_iteration_header;
    config->print_iteration = &ocp_nlp_globalization_fixed_step_print_iteration;
    config->needs_objective_value = &ocp_nlp_globalization_fixed_step_needs_objective_value;
}


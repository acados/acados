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
#include <stdlib.h>
#include <assert.h>
#include <string.h>
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_ipm.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/utils/mem.h"
// #include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_hpipm_opts_calculate_size(void *config_, void *dims_)
{
    ocp_qp_dims *dims = dims_;

    acados_size_t size = 0;
    size += sizeof(ocp_qp_hpipm_opts);
    size += sizeof(struct d_ocp_qp_ipm_arg);
    size += d_ocp_qp_ipm_arg_memsize(dims);

    size += 1 * 8;
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_hpipm_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_qp_dims *dims = dims_;
    ocp_qp_hpipm_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (ocp_qp_hpipm_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_hpipm_opts);

    opts->hpipm_opts = (struct d_ocp_qp_ipm_arg *) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_ipm_arg);

    align_char_to(8, &c_ptr);
    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_ocp_qp_ipm_arg_create(dims, opts->hpipm_opts, c_ptr);
    c_ptr += d_ocp_qp_ipm_arg_memsize(dims);

    assert((char *) raw_memory + ocp_qp_hpipm_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}


static void ocp_qp_hpipm_opts_overwrite_mode_opts(ocp_qp_hpipm_opts *opts)
{
    // overwrite some default options
    opts->hpipm_opts->res_g_max = 1e-6;
    opts->hpipm_opts->res_b_max = 1e-8;
    opts->hpipm_opts->res_d_max = 1e-8;
    opts->hpipm_opts->res_m_max = 1e-8;
    opts->hpipm_opts->iter_max = 50;
    opts->hpipm_opts->stat_max = 50;
    opts->hpipm_opts->alpha_min = 1e-8;
    opts->hpipm_opts->mu0 = 1e0;
    opts->hpipm_opts->var_init_scheme = 1;
}


void ocp_qp_hpipm_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    // ocp_qp_dims *dims = dims_;
    ocp_qp_hpipm_opts *opts = opts_;

//    d_ocp_qp_ipm_arg_set_default(SPEED, opts->hpipm_opts);
    d_ocp_qp_ipm_arg_set_default(BALANCE, opts->hpipm_opts);

    ocp_qp_hpipm_opts_overwrite_mode_opts(opts);

    return;
}



void ocp_qp_hpipm_opts_update(void *config_, void *dims_, void *opts_)
{
    //    ocp_qp_hpipm_opts *opts = (ocp_qp_hpipm_opts *)opts_;

    return;
}



void ocp_qp_hpipm_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_hpipm_opts *opts = opts_;

    const char *mode;
    if (!strcmp(field, "hpipm_mode"))
    {
        mode = (const char *) value;
        if (!strcmp(mode, "BALANCE"))
            d_ocp_qp_ipm_arg_set_default(BALANCE, opts->hpipm_opts);
        else if (!strcmp(mode, "SPEED"))
            d_ocp_qp_ipm_arg_set_default(SPEED, opts->hpipm_opts);
        else if (!strcmp(mode, "SPEED_ABS"))
            d_ocp_qp_ipm_arg_set_default(SPEED_ABS, opts->hpipm_opts);
        else if (!strcmp(mode, "ROBUST"))
            d_ocp_qp_ipm_arg_set_default(ROBUST, opts->hpipm_opts);

        ocp_qp_hpipm_opts_overwrite_mode_opts(opts);

    }
    else
    {
        d_ocp_qp_ipm_arg_set((char *) field, value, opts->hpipm_opts);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_qp_hpipm_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;
    ocp_qp_hpipm_opts *opts = opts_;

    acados_size_t size = 0;
    size += sizeof(ocp_qp_hpipm_memory);

    size += sizeof(struct d_ocp_qp_ipm_ws);

    size += d_ocp_qp_ipm_ws_memsize(dims, opts->hpipm_opts);

    size += 1 * 8;
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_qp_hpipm_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_qp_dims *dims = dims_;
    ocp_qp_hpipm_opts *opts = opts_;
    ocp_qp_hpipm_memory *mem;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (ocp_qp_hpipm_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_hpipm_memory);

    mem->hpipm_workspace = (struct d_ocp_qp_ipm_ws *) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_ipm_ws);

    struct d_ocp_qp_ipm_ws *ipm_workspace = mem->hpipm_workspace;

    align_char_to(8, &c_ptr);
    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    // ipm workspace structure
    d_ocp_qp_ipm_ws_create(dims, opts->hpipm_opts, ipm_workspace, c_ptr);
    c_ptr += ipm_workspace->memsize;

    assert((char *) raw_memory + ocp_qp_hpipm_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_qp_hpipm_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    // qp_solver_config *config = config_;
    ocp_qp_hpipm_memory *mem = mem_;

    if (!strcmp(field, "time_qp_solver_call"))
    {
        double *tmp_ptr = value;
        *tmp_ptr = mem->time_qp_solver_call;
    }
    else if (!strcmp(field, "iter"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->iter;
    }
    else if (!strcmp(field, "status"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->status;
    }
    else
    {
        printf("\nerror: ocp_qp_hpipm_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_hpipm_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    return 0;
}



/************************************************
 * functions
 ************************************************/

int ocp_qp_hpipm(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;

    qp_info *info = qp_out->misc;
    acados_timer tot_timer, qp_timer;

    acados_tic(&tot_timer);
    // cast data structures
    ocp_qp_hpipm_opts *opts = opts_;
    ocp_qp_hpipm_memory *mem = mem_;

    // zero primal solution
    // TODO add a check if warm start of first SQP iteration is implemented !!!!!!
    int ii;
    int N = qp_in->dim->N;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *ns = qp_in->dim->ns;
    for(ii=0; ii<=N; ii++)
    {
        blasfeo_dvecse(nu[ii]+nx[ii]+2*ns[ii], 0.0, qp_out->ux+ii, 0);
    }

    // solve ipm
    acados_tic(&qp_timer);
    // print_ocp_qp_in(qp_in);
    d_ocp_qp_ipm_solve(qp_in, qp_out, opts->hpipm_opts, mem->hpipm_workspace);
    d_ocp_qp_ipm_get_status(mem->hpipm_workspace, &mem->status);

    info->solve_QP_time = acados_toc(&qp_timer);
    info->interface_time = 0;  // there are no conversions for hpipm
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = mem->hpipm_workspace->iter;
    info->t_computed = 1;

    mem->time_qp_solver_call = info->solve_QP_time;
    mem->iter = mem->hpipm_workspace->iter;

    // check exit conditions
    int acados_status = mem->status;
    if (mem->status == 0) acados_status = ACADOS_SUCCESS;
    if (mem->status == 1) acados_status = ACADOS_MAXITER;
    if (mem->status == 2) acados_status = ACADOS_MINSTEP;

    return acados_status;
}


void ocp_qp_hpipm_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    // reset memory
    // void *ocp_qp_hpipm_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
    printf("acados: reset hpipm_mem\n");
    ocp_qp_hpipm_memory_assign(config_, qp_in->dim, opts_, mem_);
}

void ocp_qp_hpipm_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    ocp_qp_in *qp_in = qp_in_;
    // ocp_qp_out *qp_out = qp_out_;
    ocp_qp_hpipm_opts *opts = opts_;
    ocp_qp_hpipm_memory *mem = mem_;

    double *double_values = value;
    int nx = qp_in->dim->nx[stage];
    int nu = qp_in->dim->nu[stage];

    if (!strcmp(field, "P"))
    {
        if ((size1 != nx) || (size2 != nx))
        {
            printf("\nocp_qp_hpipm_solver_get: size of field %s not as expected, got size %d %d.\n",
                   field, size1, size2);
        }
        d_ocp_qp_ipm_get_ric_P(qp_in, opts->hpipm_opts, mem->hpipm_workspace, stage, double_values);
    }
    else if (!strcmp(field, "p"))
    {
        if ((size1 != nx) || (size2 != 1))
        {
            printf("\nocp_qp_hpipm_solver_get: size of field %s not as expected, got size %d %d.\n",
                   field, size1, size2);
        }
        d_ocp_qp_ipm_get_ric_p(qp_in, opts->hpipm_opts, mem->hpipm_workspace, stage, double_values);
    }
    else if (!strcmp(field, "K"))
    {
        if ((size1 != nu) || (size2 != nx))
        {
            printf("\nocp_qp_hpipm_solver_get: size of field %s not as expected, got size %d %d.\n",
                   field, size1, size2);
        }
        d_ocp_qp_ipm_get_ric_K(qp_in, opts->hpipm_opts, mem->hpipm_workspace, stage, double_values);
    }
    else if (!strcmp(field, "Lr"))
    {
        if ((size1 != nu) || (size2 != nu))
        {
            printf("\nocp_qp_hpipm_solver_get: size of field %s not as expected, got size %d %d.\n",
                   field, size1, size2);
        }
        d_ocp_qp_ipm_get_ric_Lr(qp_in, opts->hpipm_opts, mem->hpipm_workspace, stage, double_values);
    }
    else
    {
        printf("\nocp_qp_hpipm_solver_get: field %s not supported", field);
    }
    return;
}


void ocp_qp_hpipm_eval_sens(void *config_, void *param_qp_in_, void *sens_qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *param_qp_in = param_qp_in_;
    ocp_qp_out *sens_qp_out = sens_qp_out_;

//    qp_info *info = sens_qp_out->misc;
//    acados_timer tot_timer, qp_timer;

//    acados_tic(&tot_timer);
    // cast data structures
    ocp_qp_hpipm_opts *opts = opts_;
    ocp_qp_hpipm_memory *mem = mem_;

    // solve ipm
//    acados_tic(&qp_timer);
    // print_ocp_qp_in(param_qp_in);
    d_ocp_qp_ipm_sens(param_qp_in, sens_qp_out, opts->hpipm_opts, mem->hpipm_workspace);

//    info->solve_QP_time = acados_toc(&qp_timer);
//    info->interface_time = 0;  // there are no conversions for hpipm
//    info->total_time = acados_toc(&tot_timer);
//    info->num_iter = mem->hpipm_workspace->iter;
//    info->t_computed = 1;

    return;
}



void ocp_qp_hpipm_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = &ocp_qp_hpipm_opts_calculate_size;
    config->opts_assign = &ocp_qp_hpipm_opts_assign;
    config->opts_initialize_default = &ocp_qp_hpipm_opts_initialize_default;
    config->opts_update = &ocp_qp_hpipm_opts_update;
    config->opts_set = &ocp_qp_hpipm_opts_set;
    config->memory_calculate_size = &ocp_qp_hpipm_memory_calculate_size;
    config->memory_assign = &ocp_qp_hpipm_memory_assign;
    config->memory_get = &ocp_qp_hpipm_memory_get;
    config->workspace_calculate_size = &ocp_qp_hpipm_workspace_calculate_size;
    config->evaluate = &ocp_qp_hpipm;
    config->solver_get = &ocp_qp_hpipm_solver_get;
    config->memory_reset = &ocp_qp_hpipm_memory_reset;
    config->eval_sens = &ocp_qp_hpipm_eval_sens;

    return;
}

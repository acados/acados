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


#include "acados/ocp_nlp/ocp_nlp_cost_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_nlp_cost_dims_calculate_size(void *config_)
{
    acados_size_t size = sizeof(ocp_nlp_cost_dims);

    return size;
}



void *ocp_nlp_cost_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_dims *dims = (ocp_nlp_cost_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_dims);

    assert((char *) raw_memory + ocp_nlp_cost_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void ocp_nlp_cost_dims_set(void *config_, void *dims_, const char *field, int* value)
{
    ocp_nlp_cost_dims *dims = (ocp_nlp_cost_dims *) dims_;
    if (!strcmp(field, "nx"))
    {
        dims->nx = *value;
    }
    else if (!strcmp(field, "nz"))
    {
        dims->nz = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu = *value;
    }
    else if (!strcmp(field, "ny"))
    {
        dims->ny = *value;
    }
    else if (!strcmp(field, "ns"))
    {
        dims->ns = *value;
    }
    else if (!strcmp(field, "np"))
    {
        dims->np = *value;
    }
    else if (!strcmp(field, "np_global"))
    {
        dims->np_global = *value;
    }
    else
    {
        printf("\nerror: dimension type: %s not available in module\n", field);
        exit(1);
    }
}


void ocp_nlp_cost_dims_get(void *config_, void *dims_, const char *field, int* value)
{
    ocp_nlp_cost_dims *dims = (ocp_nlp_cost_dims *) dims_;

    if (!strcmp(field, "ny"))
    {
        *value = dims->ny;
    }
    else
    {
        printf("error: ocp_nlp_cost_dims_get: getting field %s is not implemented.\n", field);
        exit(1);
    }
}


/************************************************
 * common model
 ************************************************/

acados_size_t ocp_nlp_cost_common_model_calculate_size(ocp_nlp_cost_dims* dims)
{
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_common_model);

    size += 1 * 64;  // blasfeo_mem align
    size += 4 * blasfeo_memsize_dvec(2*ns);  // Z_usr, Z_nlp, z_usr, z_nlp

    make_int_multiple_of(8, &size);

    return size;
}



ocp_nlp_cost_common_model *ocp_nlp_cost_common_model_assign(ocp_nlp_cost_dims* dims, char **c_ptr)
{
    int ns = dims->ns;

    // struct
    ocp_nlp_cost_common_model *model = (ocp_nlp_cost_common_model *) *c_ptr;
    *c_ptr += sizeof(ocp_nlp_cost_common_model);

    // blasfeo_mem align
    align_char_to(64, c_ptr);

    // blasfeo_dvec
    assign_and_advance_blasfeo_dvec_mem(2*ns, &model->Z_usr, c_ptr);
    assign_and_advance_blasfeo_dvec_mem(2*ns, &model->Z_nlp, c_ptr);
    assign_and_advance_blasfeo_dvec_mem(2*ns, &model->z_usr, c_ptr);
    assign_and_advance_blasfeo_dvec_mem(2*ns, &model->z_nlp, c_ptr);

    // default initialization
    model->scaling = 1.0;

    return model;
}



int ocp_nlp_cost_common_model_set(ocp_nlp_cost_dims* dims, ocp_nlp_cost_common_model *model, const char *field, void *value_)
{
    int ns = dims->ns;

    if (!strcmp(field, "Zl"))
    {
        double *Zl = (double *) value_;
        blasfeo_pack_dvec(ns, Zl, 1, &model->Z_usr, 0);
    }
    else if (!strcmp(field, "Zu"))
    {
        double *Zu = (double *) value_;
        blasfeo_pack_dvec(ns, Zu, 1, &model->Z_usr, ns);
    }
    else if (!strcmp(field, "z"))
    {
        double *z = (double *) value_;
        blasfeo_pack_dvec(ns, z, 1, &model->z_usr, 0);
        blasfeo_pack_dvec(ns, z, 1, &model->z_usr, ns);
    }
    else if (!strcmp(field, "zl"))
    {
        double *zl = (double *) value_;
        blasfeo_pack_dvec(ns, zl, 1, &model->z_usr, 0);
    }
    else if (!strcmp(field, "zu"))
    {
        double *zu = (double *) value_;
        blasfeo_pack_dvec(ns, zu, 1, &model->z_usr, ns);
    }
    else if (!strcmp(field, "scaling"))
    {
        double *scaling_ptr = (double *) value_;
        model->scaling = *scaling_ptr;
    }
    else
    {
        return 0;
    }
    return 1;
}



int ocp_nlp_cost_common_model_get(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, const char *field, void *value_)
{
    int ns = dims->ns;

    double *value = (double *) value_;

    if (!strcmp(field, "Zl"))
    {
        blasfeo_unpack_dvec(ns, &model->Z_usr, 0, value, 1);
    }
    else if (!strcmp(field, "Zu"))
    {
        blasfeo_unpack_dvec(ns, &model->Z_usr, ns, value, 1);
    }
    else if (!strcmp(field, "zl"))
    {
        blasfeo_unpack_dvec(ns, &model->z_usr, 0, value, 1);
    }
    else if (!strcmp(field, "zu"))
    {
        blasfeo_unpack_dvec(ns, &model->z_usr, ns, value, 1);
    }
    else if (!strcmp(field, "scaling"))
    {
        value[0] = model->scaling;
    }
    else
    {
        return 0;
    }
    return 1;
}

void cost_common_add_slack_contributions_to_fun_and_scale(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, ocp_nlp_cost_common_memory *memory, struct blasfeo_dvec *tmp_2ns)
{
    int ns = dims->ns;

    struct blasfeo_dvec *ux = memory->ux;

    // slack update function value
    // tmp_2ns = 2 * z + Z .* slack
    blasfeo_dveccpsc(2*ns, 2.0, &model->z_nlp, 0, tmp_2ns, 0);
    blasfeo_dvecmulacc(2*ns, &model->Z_nlp, 0, ux, dims->nx+dims->nu, tmp_2ns, 0);
    // fun += .5 * (tmp_2ns .* slack)
    memory->fun += 0.5 * blasfeo_ddot(2*ns, tmp_2ns, 0, ux, dims->nx+dims->nu);

    // scale
    memory->fun *= model->scaling;
}

void cost_common_update_gradient_with_slacks_and_scale(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, ocp_nlp_cost_common_memory *memory)
{
    int ns = dims->ns;
    int nx = dims->nx;
    int nu = dims->nu;
    struct blasfeo_dvec *ux = memory->ux;

    // slack update gradient
    blasfeo_dveccp(2*ns, &model->z_nlp, 0, &memory->grad, nx+nu);
    blasfeo_dvecmulacc(2*ns, &model->Z_nlp, 0, ux, nx+nu, &memory->grad, nx+nu);

    // scale
    if (model->scaling != 1.0)
    {
        blasfeo_dvecsc(nu+nx+2*ns, model->scaling, &memory->grad, 0);
    }
}

void ocp_nlp_cost_common_initialize(ocp_nlp_cost_dims *dims, ocp_nlp_cost_common_model *model, ocp_nlp_cost_common_memory *memory)
{
    int ns = dims->ns;

    // adjust according to orphan_mask
    if (memory->orphan_mask)
    {
        // z_nlp
        blasfeo_dvecmul(2*ns, &model->z_usr, 0, memory->orphan_mask, 0, &model->z_nlp, 0);

        // Z_nlp
        for (int ii = 0; ii < 2*ns; ii++)
            if (BLASFEO_DVECEL(memory->orphan_mask, ii) == 0)
                BLASFEO_DVECEL(&model->Z_nlp, ii) = 1.0;
            else
                BLASFEO_DVECEL(&model->Z_nlp, ii) = BLASFEO_DVECEL(&model->Z_usr, ii);
    }
    else
    {
        // z: just pointer alias
        model->Z_nlp = model->Z_usr;
        model->z_nlp = model->z_usr;
    }

    // Z_qp
    blasfeo_dveccpsc(2*ns, model->scaling, &model->Z_nlp, 0, memory->Z, 0);

    return;
}


/************************************************
 * common memory
 ************************************************/

acados_size_t ocp_nlp_cost_common_memory_calculate_size(ocp_nlp_cost_dims *dims)
{
    int ns = dims->ns;
    int nx = dims->nx;
    int nu = dims->nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_common_memory);

    size += sizeof(ocp_nlp_cost_capsule);  // cost capsule

    size += 1 * 64;  // blasfeo_mem align
    size += 1 * blasfeo_memsize_dvec(nu + nx + 2*ns);  // grad

    make_int_multiple_of(8, &size);

    return size;
}



ocp_nlp_cost_common_memory *ocp_nlp_cost_common_memory_assign(ocp_nlp_cost_dims *dims, char **c_ptr)
{
    int ns = dims->ns;
    int nx = dims->nx;
    int nu = dims->nu;

    // struct
    ocp_nlp_cost_common_memory *memory = (ocp_nlp_cost_common_memory *) *c_ptr;
    *c_ptr += sizeof(ocp_nlp_cost_common_memory);

    // cost capsule (bundling of module pointers, filled in precompute)
    memory->capsule = (ocp_nlp_cost_capsule *) *c_ptr;
    *c_ptr += sizeof(ocp_nlp_cost_capsule);

    // blasfeo_mem align
    align_char_to(64, c_ptr);

    // blasfeo_dvec
    assign_and_advance_blasfeo_dvec_mem(nu + nx + 2*ns, &memory->grad, c_ptr);

    memory->orphan_mask = NULL;

    return memory;
}


int ocp_nlp_cost_common_memory_set(ocp_nlp_cost_common_memory *memory, const char *field, void *value)
{
    if (!strcmp(field, "ux_ptr"))
    {
        memory->ux = value;
    }
    else if (!strcmp(field, "z_alg_ptr"))
    {
        memory->z_alg = value;
    }
    else if (!strcmp(field, "dzdux_tran_ptr"))
    {
        memory->dzdux_tran = value;
    }
    else if (!strcmp(field, "RSQrq_ptr"))
    {
        memory->RSQrq = value;
    }
    else if (!strcmp(field, "Z_ptr"))
    {
        memory->Z = value;
    }
    else if (!strcmp(field, "orphan_mask_ptr"))
    {
        memory->orphan_mask = value;
    }
    else if (!strcmp(field, "jac_lag_stat_p_global_ptr"))
    {
        memory->jac_lag_stat_p_global = value;
    }
    else if (!strcmp(field, "adj_lag_p_global_ptr"))
    {
        memory->adj_lag_p_global = value;
    }
    else if (!strcmp(field, "seed_ux_ptr"))
    {
        memory->seed_ux = value;
    }
    else
    {
        return 0;
    }
    return 1;
}


// returns pointer to the requested field if it is handled by the common memory,
// NULL if the field is not available there (e.g. module-specific fields).
void *ocp_nlp_cost_common_memory_get(ocp_nlp_cost_common_memory *memory, const char *field)
{
    if (!strcmp(field, "fun"))
    {
        return &memory->fun;
    }
    else if (!strcmp(field, "grad"))
    {
        return &memory->grad;
    }
    else if (!strcmp(field, "cost_capsule"))
    {
        return memory->capsule;
    }
    else
    {
        return NULL;
    }
}


// fills the cost capsule in the common cost memory with the module pointers.
// Intended to be called from the precompute function of the cost modules,
// where all module pointers are available.
void ocp_nlp_cost_common_fill_capsule(ocp_nlp_cost_common_memory *common_memory, void *config_, void *dims_, void *model_, void *opts_,
                                      void *memory_, void *work_)
{
    common_memory->capsule->config = config_;
    common_memory->capsule->dims = dims_;
    common_memory->capsule->model = model_;
    common_memory->capsule->opts = opts_;
    common_memory->capsule->memory = memory_;
    common_memory->capsule->work = work_;

    return;
}



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_cost_common_opts_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_common_opts);
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_cost_common_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_common_opts *opts = (ocp_nlp_cost_common_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_common_opts);

    assert((char *) raw_memory + ocp_nlp_cost_common_opts_calculate_size(config_, dims_) >= c_ptr);

    return opts;
}



void ocp_nlp_cost_common_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_cost_common_opts *opts = opts_;

    opts->compute_hess = 1;
    opts->exact_hess = 0;
    opts->use_numerical_hessian = 0;
    opts->integrator_cost = 0;
    opts->with_solution_sens_wrt_params_forw = 0;
    opts->with_solution_sens_wrt_params_adj = 0;
    opts->add_hess_contribution = 0;

    return;
}



void ocp_nlp_cost_common_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_nlp_cost_common_opts *opts = opts_;

    if (!strcmp(field, "exact_hess"))
    {
        int *int_ptr = value;
        opts->exact_hess = *int_ptr;
    }
    else if (!strcmp(field, "compute_hess"))
    {
        int *int_ptr = value;
        opts->compute_hess = *int_ptr;
    }
    else if (!strcmp(field, "numerical_hessian"))
    {
        int *opt_val = (int *) value;
        opts->use_numerical_hessian = *opt_val;
    }
    else if (!strcmp(field, "integrator_cost"))
    {
        int *opt_val = (int *) value;
        opts->integrator_cost = *opt_val;
    }
    else if (!strcmp(field, "with_solution_sens_wrt_params_forw"))
    {
        int *opt_val = (int *) value;
        opts->with_solution_sens_wrt_params_forw = *opt_val;
    }
    else if (!strcmp(field, "with_solution_sens_wrt_params_adj"))
    {
        int *opt_val = (int *) value;
        opts->with_solution_sens_wrt_params_adj = *opt_val;
    }
    else if (!strcmp(field, "add_hess_contribution"))
    {
        int *int_ptr = value;
        opts->add_hess_contribution = *int_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_cost_common_opts_set\n", field);
        exit(1);
    }

    return;
}



int ocp_nlp_cost_common_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    ocp_nlp_cost_common_opts *opts = opts_;

    if (!strcmp(field, "add_hess_contribution"))
    {
        int *int_ptr = value;
        *int_ptr = opts->add_hess_contribution;
    }
    else
    {
        return 0;
    }
    return 1;
}



/************************************************
 * config
 ************************************************/

acados_size_t ocp_nlp_cost_config_calculate_size()
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_config);

    return size;
}

ocp_nlp_cost_config *ocp_nlp_cost_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_nlp_cost_config *config = (ocp_nlp_cost_config *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_config);

    return config;
}

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

// This is a template based custom_update function

#include <stdlib.h>
#include <stdio.h>

#include "acados_solver_{{ name }}.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/mem.h"

// #include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

#include "helpers_{{ name }}.h"
#include "{{ name }}_model/{{ name }}_model.h"

typedef struct custom_memory
{
    external_function_casadi* ext_helpers;
    void *raw_memory; // Pointer to allocated memory, to be used for freeing
} custom_memory;


static void *example_custom_memory_create({{ name }}_solver_capsule* capsule)
{
    static external_function_casadi ext_helpers;
    ext_helpers.casadi_fun = &helpers;
    ext_helpers.casadi_work = &helpers_work;
    ext_helpers.casadi_sparsity_in = &helpers_sparsity_in;
    ext_helpers.casadi_sparsity_out = &helpers_sparsity_out;
    ext_helpers.casadi_n_in = &helpers_n_in;
    ext_helpers.casadi_n_out = &helpers_n_out;

    acados_size_t bytes = sizeof(custom_memory);
    bytes += external_function_casadi_calculate_size(&ext_helpers);

    void *ptr = acados_calloc(1, bytes);
    char *c_ptr = (char *) ptr;
    custom_memory *custom_mem = (custom_memory *) c_ptr;
    custom_mem->ext_helpers = &ext_helpers;

    c_ptr += sizeof(custom_memory);
    external_function_casadi_assign(custom_mem->ext_helpers, c_ptr);
    c_ptr += external_function_casadi_calculate_size(&ext_helpers);

    custom_mem->raw_memory = ptr;

    return custom_mem;
}


int custom_update_init_function({{ name }}_solver_capsule* capsule)
{
    capsule->custom_update_memory = example_custom_memory_create(capsule);
    return 1;
}


int custom_update_function({{ name }}_solver_capsule* capsule, double* data, int data_len)
{
    custom_memory *custom_mem = (custom_memory *) capsule->custom_update_memory;
    external_function_casadi* fun = custom_mem->ext_helpers;
    fun->args[0] = data;
    // TODO! Ensure data_len == len(p_slow)

{% set n_pools = casadi_pool_names | length %}
{% for ip in range(end=n_pools) %}
{% set pool_name = casadi_pool_names[ip] %}
{% set fun_name_split = pool_name | split(pat='|') %}
    fun->res[{{ ip }}] = {{ fun_name_split[0] }}_get_pool_double("{{ pool_name }}");
{%- endfor %}

    fun->casadi_fun((const double **) fun->args, fun->res, fun->iw, fun->w, NULL);
    return 1;
}

int custom_update_terminate_function({{ name }}_solver_capsule* capsule)
{
    custom_memory *mem = capsule->custom_update_memory;
    free(mem->raw_memory);
    return 1;
}

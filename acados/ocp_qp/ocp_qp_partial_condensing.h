/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_H_
#define ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp_red.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"



typedef struct
{
    ocp_qp_dims *orig_dims;
    ocp_qp_dims *red_dims; // dims of reduced qp
    ocp_qp_dims *pcond_dims;
    int *block_size;
    int N2;
    int N2_bkp;
} ocp_qp_partial_condensing_dims;



typedef struct ocp_qp_partial_condensing_opts_
{
    struct d_part_cond_qp_arg *hpipm_pcond_opts;
    struct d_ocp_qp_reduce_eq_dof_arg *hpipm_red_opts;
    int N2;
    int N2_bkp;
    int *block_size;
    bool block_size_was_set;
    int ric_alg;
    int mem_qp_in; // allocate qp_in in memory
} ocp_qp_partial_condensing_opts;



typedef struct ocp_qp_partial_condensing_memory_
{
    struct d_part_cond_qp_ws *hpipm_pcond_work;
    struct d_ocp_qp_reduce_eq_dof_ws *hpipm_red_work;
    // in memory
    ocp_qp_in *pcond_qp_in;
    ocp_qp_out *pcond_qp_out;
    ocp_qp_seed *pcond_qp_seed;
    ocp_qp_in *red_qp; // reduced qp
    ocp_qp_out *red_sol; // reduced qp sol
    ocp_qp_seed *red_seed;
    // only pointer
    ocp_qp_in *ptr_qp_in;
    ocp_qp_in *ptr_pcond_qp_in;
    ocp_qp_seed *ptr_qp_seed;
    qp_info *qp_out_info; // info in pcond_qp_in
    ocp_qp_partial_condensing_dims *dims;
    double time_qp_xcond;
} ocp_qp_partial_condensing_memory;



//
acados_size_t ocp_qp_partial_condensing_opts_calculate_size(void *dims);
//
void *ocp_qp_partial_condensing_opts_assign(void *dims, void *raw_memory);
//
void ocp_qp_partial_condensing_opts_initialize_default(void *dims, void *opts_);
//
void ocp_qp_partial_condensing_opts_update(void *dims, void *opts_);
//
void ocp_qp_partial_condensing_opts_set(void *opts_, const char *field, void* value);
//
acados_size_t ocp_qp_partial_condensing_memory_calculate_size(void *dims, void *opts_);
//
void *ocp_qp_partial_condensing_memory_assign(void *dims, void *opts, void *raw_memory);
//
acados_size_t ocp_qp_partial_condensing_workspace_calculate_size(void *dims, void *opts_);
//
int ocp_qp_partial_condensing(void *in, void *out, void *opts, void *mem, void *work);
//
int ocp_qp_partial_condensing_condense_lhs(void *in, void *out, void *opts, void *mem, void *work);
//
int ocp_qp_partial_condensing_condense_rhs(void *in, void *out, void *opts, void *mem, void *work);
//
int ocp_qp_partial_expansion(void *in, void *out, void *opts, void *mem, void *work);
//
void ocp_qp_partial_condensing_config_initialize_default(void *config_);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_H_

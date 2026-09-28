/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_SOLVER_H_
#define ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"



typedef struct
{
    ocp_qp_dims *orig_dims;
    void *xcond_dims;
} ocp_qp_xcond_solver_dims;



typedef struct ocp_qp_xcond_solver_opts_
{
    void *xcond_opts;
    void *qp_solver_opts;
    bool initialize_next_xcond_qp_from_qp_out;
} ocp_qp_xcond_solver_opts;



typedef struct ocp_qp_xcond_solver_memory_
{
    void *xcond_memory;
    void *solver_memory;
    void *xcond_qp_in;
    void *xcond_qp_out;
    void *xcond_seed;
} ocp_qp_xcond_solver_memory;



typedef struct ocp_qp_xcond_solver_workspace_
{
    void *xcond_work;
    void *qp_solver_work;
} ocp_qp_xcond_solver_workspace;



typedef struct
{
    acados_size_t (*dims_calculate_size)(void *config, int N);
    ocp_qp_xcond_solver_dims *(*dims_assign)(void *config, int N, void *raw_memory);
    void (*dims_set)(void *config_, ocp_qp_xcond_solver_dims *dims, int stage, const char *field, int* value);
    void (*dims_get)(void *config_, ocp_qp_xcond_solver_dims *dims, int stage, const char *field, int* value);
    acados_size_t (*opts_calculate_size)(void *config, ocp_qp_xcond_solver_dims *dims);
    void *(*opts_assign)(void *config, ocp_qp_xcond_solver_dims *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, ocp_qp_xcond_solver_dims *dims, void *opts);
    void (*opts_update)(void *config, ocp_qp_xcond_solver_dims *dims, void *opts);
    void (*opts_set)(void *config_, void *opts_, const char *field, void* value);
    void (*opts_get)(void *config_, void *opts_, const char *field, void* value);
    acados_size_t (*memory_calculate_size)(void *config, ocp_qp_xcond_solver_dims *dims, void *opts);
    void *(*memory_assign)(void *config, ocp_qp_xcond_solver_dims *dims, void *opts, void *raw_memory);
    void (*memory_get)(void *config_, void *mem_, const char *field, void* value);
    void (*solver_get)(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
    void (*memory_reset)(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts, void *mem, void *work);
    acados_size_t (*workspace_calculate_size)(void *config, ocp_qp_xcond_solver_dims *dims, void *opts);
    int (*evaluate)(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts, void *mem, void *work);
    int (*condense_lhs)(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts, void *mem, void *work);
    int (*condense_rhs_and_solve)(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts, void *mem, void *work);
    void (*eval_forw_sens)(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_seed *seed, ocp_qp_out *sens_qp_out, void *opts, void *mem, void *work);
    void (*eval_adj_sens)(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_seed *seed, ocp_qp_out *sens_qp_out, void *opts, void *mem, void *work);
    void (*terminate)(void *config, void *mem, void *work);
    qp_solver_config *qp_solver;  // either ocp_qp_solver or dense_solver
    ocp_qp_xcond_config *xcond;
} ocp_qp_xcond_solver_config;  // pcond - partial condensing or fcond - full condensing


typedef struct ocp_qp_xcond_solver
{
    ocp_qp_xcond_solver_config *config;
    ocp_qp_xcond_solver_dims *dims;
    ocp_qp_xcond_solver_opts *opts;
    ocp_qp_xcond_solver_memory *mem;
    ocp_qp_xcond_solver_workspace *work;
} ocp_qp_xcond_solver;



/* config */
//
acados_size_t ocp_qp_xcond_solver_config_calculate_size();
//
ocp_qp_xcond_solver_config *ocp_qp_xcond_solver_config_assign(void *raw_memory);

/* dims */
//
acados_size_t ocp_qp_xcond_solver_dims_calculate_size(void *config, int N);
//
ocp_qp_xcond_solver_dims *ocp_qp_xcond_solver_dims_assign(void *config, int N, void *raw_memory);
//
void ocp_qp_xcond_solver_dims_set_(void *config, ocp_qp_xcond_solver_dims *dims, int stage, const char *field, int* value);

/* opts */
//
acados_size_t ocp_qp_xcond_solver_opts_calculate_size(void *config, ocp_qp_xcond_solver_dims *dims);
//
void *ocp_qp_xcond_solver_opts_assign(void *config, ocp_qp_xcond_solver_dims *dims, void *raw_memory);
//
void ocp_qp_xcond_solver_opts_initialize_default(void *config, ocp_qp_xcond_solver_dims *dims, void *opts_);
//
void ocp_qp_xcond_solver_opts_update(void *config, ocp_qp_xcond_solver_dims *dims, void *opts_);
//
void ocp_qp_xcond_solver_opts_set_(void *config_, void *opts_, const char *field, void* value);
void ocp_qp_xcond_solver_opts_get_(void *config_, void *opts_, const char *field, void* value);

/* memory */
//
acados_size_t ocp_qp_xcond_solver_memory_calculate_size(void *config, ocp_qp_xcond_solver_dims *dims, void *opts_);
//
void *ocp_qp_xcond_solver_memory_assign(void *config, ocp_qp_xcond_solver_dims *dims, void *opts_, void *raw_memory);

/* workspace */
//
acados_size_t ocp_qp_xcond_solver_workspace_calculate_size(void *config, ocp_qp_xcond_solver_dims *dims, void *opts_);

/* config */
//
int ocp_qp_xcond_solve(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *mem_, void *work_);

int ocp_qp_xcond_cond_lhs(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *mem_, void *work_);

int ocp_qp_xcond_cond_rhs_and_solve(void *config, ocp_qp_xcond_solver_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *mem_, void *work_);


//
void ocp_qp_xcond_solver_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_PARTIAL_CONDENSING_SOLVER_H_

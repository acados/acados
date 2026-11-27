/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_OCP_QP_OCP_QP_DAQP_H_
#define ACADOS_OCP_QP_OCP_QP_DAQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// daqp
#include "daqp/include/types.h"
#include "daqp/include/multistage.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"


/*
 * Structure-exploiting DAQP (daqp/include/multistage.h, with its BLASFEO backend) for
 * the OCP QP, without condensing. DAQP uses the matrices of qp_in (BAbt,
 * RSQrq, DCt) directly; the vectors are copied. What has changed since the
 * previous call (compared with copies) decides how much of the factorization
 * DAQP updates, and the working set is kept as a warm start.
 *
 * Soft constraints are soft constraints of DAQP (each slack has to soften a
 * single constraint): a slack s >= s_min with the penalty 0.5*Z*s^2 + z*s is
 * the nonnegative slack t = s - s_hat of DAQP, with s_hat = max(s_min, -z/Z),
 * the bound shifted by s_hat, and the weights rho = 1/Z and w = z + Z*s_hat.
 *
 * Binary constraints (solved by DAQP's branch and bound): binary is 0/1 per
 * constraint, in the order of DAQP's sense. binary_horizon: -1 (default) keeps
 * them at all stages, k > 0 only at the stages t < k, and 0 relaxes them to
 * their ranges (with multipliers; switching to and from 0 needs no new setup).
 * After branch and bound, lam and pi are zero.
 * rel_subopt and abs_subopt (DAQP's settings, default 0) stop branch and bound
 * early: a node has to improve the objective J of the best binary solution by
 * more than abs_subopt + rel_subopt*|J|.
 */

typedef struct ocp_qp_daqp_opts_
{
    DAQPSettings *daqp_opts;
    int warm_start;
    int print_level;
    int *binary;  // 0/1 per constraint
    int m_binary;  // length of binary
    int binary_horizon;  // -1: all stages, 0: relaxed, k > 0: stages t < k
} ocp_qp_daqp_opts;



typedef struct ocp_qp_daqp_memory_
{
    DAQPMultistage ocp;  // the QP in the format of DAQP (pointers to the data below)
    DAQPWorkspace *work;
    void *daqp_mem;  // memory of DAQP's multistage solver (worst case of the dimensions)
    acados_size_t daqp_mem_size;
    int setup_done;  // work is set up for ocp

    double *z;  // [x_0; u_0; x_1; ...] (format of DAQP)
    double *lam;  // multipliers of the bounds and general constraints (signed)
    double *pi;  // multipliers of the dynamics
    // copies of the matrices of the previous call (the first nu+nx rows of
    // BAbt and RSQrq, and DCt), to detect changes
    struct blasfeo_dmat *cBAbt, *cRSQrq, *cDCt;

    // soft constraints (per constraint)
    int *sense;
    double *rho_l, *rho_u, *w_l, *w_u;  // weights of DAQP
    double *s_l, *s_u;  // shifts s_hat of the bounds

    double *xold;  // for branch and bound
    double *u;  // the array u of the workspace (branch and bound swaps work->u and work->xold)

    double time_qp_solver_call;
    int iter;
    int status;
} ocp_qp_daqp_memory;



acados_size_t ocp_qp_daqp_opts_calculate_size(void *config, void *dims);
//
void *ocp_qp_daqp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_qp_daqp_opts_initialize_default(void *config, void *dims, void *opts_);
//
void ocp_qp_daqp_opts_update(void *config, void *dims, void *opts_);
//
void ocp_qp_daqp_opts_set(void *config_, void *opts_, const char *field, void *value);
//
void ocp_qp_daqp_opts_get(void *config_, void *opts_, const char *field, void *value);
//
acados_size_t ocp_qp_daqp_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_qp_daqp_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
void ocp_qp_daqp_memory_get(void *config_, void *mem_, const char *field, void* value);
//
acados_size_t ocp_qp_daqp_workspace_calculate_size(void *config, void *dims, void *opts_);
//
int ocp_qp_daqp(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_daqp_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_);
//
void ocp_qp_daqp_eval_forw_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_daqp_eval_adj_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_daqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2);
//
void ocp_qp_daqp_terminate(void *config_, void *mem_, void *work_);
//
void ocp_qp_daqp_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_DAQP_H_

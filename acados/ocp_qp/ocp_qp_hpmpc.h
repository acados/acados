/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_OCP_QP_OCP_QP_HPMPC_H_
#define ACADOS_OCP_QP_OCP_QP_HPMPC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef enum hpmpc_options_t_ { HPMPC_DEFAULT_ARGUMENTS } hpmpc_options_t;

typedef struct ocp_qp_hpmpc_opts_
{
    double tol;
    int max_iter;
    double mu0;
    double alpha_min;
    int warm_start;
    int N2;  // horizion length of the partially condensed problem

    // partial tightening
    double sigma_mu;
    int N;
    int M;
} ocp_qp_hpmpc_opts;

// struct of the solver memory
typedef struct ocp_qp_hpmpc_memory_
{
    struct blasfeo_dvec *hpi;
    double *stats;

    // workspace
    void *hpmpc_work;  // raw workspace

    // partial tightening-specific (init of extra variables)
    struct blasfeo_dvec *lam0;
    struct blasfeo_dvec *ux0;
    struct blasfeo_dvec *pi0;
    struct blasfeo_dvec *t0;

    // 2. workspace
    struct blasfeo_dmat *hsL;
    struct blasfeo_dmat *hsric_work_mat;
    struct blasfeo_dmat sLxM;
    struct blasfeo_dmat sPpM;

    struct blasfeo_dvec *hsQx;
    struct blasfeo_dvec *hsqx;
    struct blasfeo_dvec *hstinv;
    struct blasfeo_dvec *hsrq;
    struct blasfeo_dvec *hsdux;

    struct blasfeo_dvec *hsdlam;
    struct blasfeo_dvec *hsdt;
    struct blasfeo_dvec *hsdpi;
    struct blasfeo_dvec *hslamt;

    struct blasfeo_dvec *hsPb;

    void *work_ric;

    int out_iter;

    double time_qp_solver_call;
    int iter;
    int status;

} ocp_qp_hpmpc_memory;

acados_size_t ocp_qp_hpmpc_opts_calculate_size(void *config_, ocp_qp_dims *dims);
//
void *ocp_qp_hpmpc_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_hpmpc_opts_initialize_default(void *config_, ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_hpmpc_opts_update(void *config_, ocp_qp_dims *dims, void *opts_);
//
acados_size_t ocp_qp_hpmpc_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_);
//
void *ocp_qp_hpmpc_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_, void *raw_memory);
//
acados_size_t ocp_qp_hpmpc_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_);
//
int ocp_qp_hpmpc(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_hpmpc_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_);
//
void ocp_qp_hpmpc_eval_forw_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_hpmpc_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
//

//
void ocp_qp_hpmpc_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_HPMPC_H_

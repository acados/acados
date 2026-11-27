#ifndef ACADOS_DENSE_LP_FEASIBILITY_DENSE_LP_FEASIBILITY_H_
#define ACADOS_DENSE_LP_FEASIBILITY_DENSE_LP_FEASIBILITY_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

/* Forward declaration to avoid pulling simplex internals into the public header */
struct d_simplex_workspace;

#define DENSE_LP_PHASE1_ARTIFICIAL 0
#define DENSE_LP_PHASE1_NAZARETH 1

typedef struct dense_lp_feasibility_opts_
{
    int print_level; // 0: no printing, 1: print warning & final report, 2: print LP in, print simplex iter and basis
    int max_iter;
    int phase1_method;
    double tol;
    double tol_infeas;
    double tol_pivot;
} dense_lp_feasibility_opts;

typedef struct dense_lp_feasibility_memory_
{
    struct d_simplex_workspace *simplex_workspace;
    double time_lp_solver_call;  // setup + solve
    double time_lp_setup;        // building the Phase-I LP from the dense QP
    int iter;
    int status;
} dense_lp_feasibility_memory;

// opts
acados_size_t dense_lp_feasibility_opts_calculate_size(void *config, void *dims);
void *dense_lp_feasibility_opts_assign(void *config, void *dims, void *raw_memory);
void dense_lp_feasibility_opts_initialize_default(void *config, void *dims, void *opts_);
void dense_lp_feasibility_opts_update(void *config, void *dims, void *opts_);
void dense_lp_feasibility_opts_set(void *config, void *opts_, const char *field, void *value);
void dense_lp_feasibility_opts_get(void *config, void *opts_, const char *field, void *value);

// memory
acados_size_t dense_lp_feasibility_memory_calculate_size(void *config_, void *dims_, void *opts_);
void *dense_lp_feasibility_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory);
void dense_lp_feasibility_memory_get(void *config_, void *mem_, const char *field, void *value);

// workspace
acados_size_t dense_lp_feasibility_workspace_calculate_size(void *config_, void *dims_, void *opts_);

// solver
int dense_lp_feasibility(void *config, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_);

// solver config
void dense_lp_feasibility_config_initialize_default(void *config_);

// misc
void dense_lp_feasibility_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work);
void dense_lp_feasibility_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_,
                                     const char *field, int stage, void *value, int size1, int size2);
void dense_lp_feasibility_terminate(void *config_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_DENSE_LP_FEASIBILITY_DENSE_LP_FEASIBILITY_H_
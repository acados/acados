/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef ACADOS_OCP_NLP_OCP_NLP_GN_SQP_H_
#define ACADOS_OCP_NLP_OCP_NLP_GN_SQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"

#ifdef YT
#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_rk_common_yt.h"
#else
#include "acados/sim/sim_collocation.h"
#include "acados/sim/sim_rk_common.h"
#endif

#include "acados/utils/types.h"
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"

typedef struct {
    // general nlp options
    int maxIter;
    // specific SQP options
    ocp_qp_xcond_solver *qp_solver;
    #ifdef YT
    sim_solver_yt **sim_solvers;
    #else
    #endif
    // QP solver options
    void *qp_solver_args;
    #ifdef YT
    // integrator options
    void **sim_solvers_args;
    #else
    #endif
} ocp_nlp_gn_sqp_args;



typedef struct {
    int num_vars;
    double **x;
    double **u;
    double **pi;
    double **lam;

    ocp_nlp_dims *dims;
    void *qp_solver_mem;
} ocp_nlp_gn_sqp_memory;



typedef struct {
    double *w;

    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;
    void *qp_work;

    struct d_strvec *tmp_vecs;  // N+1 vectors of dimension nx[i]+nu[i] to store interm. results
                                // not using max(nx+nu) for parallelization in the future

} ocp_nlp_gn_sqp_work;


//
#ifdef YT
int ocp_nlp_gn_sqp_calculate_args_size(ocp_nlp_dims *dims, qp_solver_t qp_solver_name, sim_solver_t *sim_solver_names, int *num_stages);
#else
int ocp_nlp_gn_sqp_calculate_args_size(ocp_nlp_dims *dims, qp_solver_t qp_solver_name);
#endif
//
#ifdef YT
ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_assign_args(ocp_nlp_dims *dims, qp_solver_t qp_solver_name,
     sim_solver_t *sim_solver_names, int *num_stages, void *raw_memory);
#else
ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_assign_args(ocp_nlp_dims *dims, qp_solver_t qp_solver_name, void *raw_memory);
#endif
//
int ocp_nlp_gn_sqp_calculate_memory_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args);
//
ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_assign_memory(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args, void *raw_memory);
//
int ocp_nlp_gn_sqp_calculate_workspace_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args);
//
int ocp_nlp_gn_sqp(ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, void *work_);
//

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_GN_SQP_H_

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

// acados_c
#include <acados_c/ocp_qp.h>
#include <acados_c/sim.h>
// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_collocation.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"

typedef struct
{
    int maxIter;
	double min_res_g;
	double min_res_b;
	double min_res_d;
	double min_res_m;
    // QP solver
    ocp_qp_xcond_solver_fcn_ptrs *qp_solver;
    void *qp_solver_args;
    // integrators
    sim_solver_fcn_ptrs **sim_solvers;
    void **sim_solvers_args;
} ocp_nlp_gn_sqp_args;



typedef struct
{
    ocp_nlp_dims *dims;
    void *qp_solver_mem;

    void **sim_solvers_mem;

    // residuals
	ocp_nlp_res *nlp_res;

	// nlp memory
	ocp_nlp_mem *nlp_mem;

	int sqp_iter;

} ocp_nlp_gn_sqp_memory;



typedef struct
{

    // QP solver
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;
    void *qp_work;

    // integrators
    sim_in **sim_in;
    sim_out **sim_out;
    void **sim_solvers_work;

	// temporary stuff
    // N+1 vectors of dimension nx[i]+nu[i] to store interm. results
    // not using max(nx+nu) for parallelization in the future
	// XXX take Max instead ?????
	struct blasfeo_dmat *tmp_ny_ny;
	struct blasfeo_dmat *tmp_nv_ny;
	struct blasfeo_dvec *tmp_nbg;
    struct blasfeo_dvec *tmp_nux;

} ocp_nlp_gn_sqp_work;



int ocp_nlp_gn_sqp_calculate_args_size(ocp_nlp_dims *dims, ocp_qp_xcond_solver_fcn_ptrs *qp_solver, sim_solver_fcn_ptrs *sim_solvers);
//
ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_assign_args(ocp_nlp_dims *dims, ocp_qp_xcond_solver_fcn_ptrs *qp_solver, sim_solver_fcn_ptrs *sim_solvers, void *raw_memory);
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

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

#ifndef ACADOS_OCP_QP_OCP_QP_COMMON_H_
#define ACADOS_OCP_QP_OCP_QP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
#include "hpipm/include/hpipm_d_ocp_qp_res.h"
// acados
#include "acados/utils/types.h"


typedef struct d_ocp_qp_dim ocp_qp_dims;
typedef struct d_ocp_qp ocp_qp_in;
typedef struct d_ocp_qp_sol ocp_qp_out;
typedef struct d_ocp_qp_res ocp_qp_res;
typedef struct d_ocp_qp_res_workspace ocp_qp_res_ws;



typedef struct
{
	// TODO
} ocp_qp_dims_stage;



typedef struct
{
	struct blasfeo_dmat *BAbt;
	struct blasfeo_dvec *b;
	struct blasfeo_dmat *RSQrq;
	struct blasfeo_dvec *rq;
	struct blasfeo_dmat *DCt;
	struct blasfeo_dvec *d;
	struct blasfeo_dvec *Z;
	struct blasfeo_dvec *z;
	int **idxb;
	int **idxs;
} ocp_qp_in_stage;



#ifndef QP_SOLVER_CONFIG_
#define QP_SOLVER_CONFIG_

typedef struct
{
    int (*evaluate) (void *config, void *qp_in, void *qp_out, void *args, void *mem, void *work);
    int (*opts_calculate_size) (void *config, void *dims);
    void *(*opts_assign) (void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *args);
    int (*memory_calculate_size)(void *config, void *dims, void *args);
    void *(*memory_assign)(void *config, void *dims, void *args, void *raw_memory);
    int (*workspace_calculate_size)(void *config, void *dims, void *args);
} qp_solver_config;

#endif



typedef struct
{
    int (*evaluate) (void *config, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem, void *work);
    int (*opts_calculate_size) (void *config, ocp_qp_dims *dims);
    void *(*opts_assign) (void *config, ocp_qp_dims *dims, void *raw_memory);
    void *(*copy_args) (ocp_qp_dims *dims, void *raw_memory, void *source_); // ???
    void (*opts_initialize_default) (void *config, void *args);
    int (*memory_calculate_size) (void *config, ocp_qp_dims *dims, void *args);
    void *(*memory_assign) (void *config, ocp_qp_dims *dims, void *args, void *raw_memory);
    int (*workspace_calculate_size) (void *config, ocp_qp_dims *dims, void *args);
    qp_solver_config *qp_solver; // either ocp_qp_solver or dense_solver
	int N2;
} ocp_qp_xcond_solver_config;




typedef struct
{
    double solve_QP_time;
    double condensing_time;
    double interface_time;
    double total_time;
    int    num_iter;
} ocp_qp_info;


//
int ocp_qp_solver_config_calculate_size();
//
qp_solver_config *ocp_qp_solver_config_assign(void *raw_memory);
//
int ocp_qp_xcond_solver_config_calculate_size();
//
ocp_qp_xcond_solver_config *ocp_qp_xcond_solver_config_assign(void *raw_memory);
//
int ocp_qp_dims_calculate_size(int N);
//
ocp_qp_dims *ocp_qp_dims_assign(int N, void *raw_memory);
//
int ocp_qp_dims_stage_calculate_size();
//
ocp_qp_dims_stage *ocp_qp_dims_stage_assign(void *raw_memory);
//
int ocp_qp_in_calculate_size(void *config, ocp_qp_dims *dims);
//
ocp_qp_in *ocp_qp_in_assign(void *config, ocp_qp_dims *dims, void *raw_memory);
//
int ocp_qp_in_stage_calculate_size(void *config, ocp_qp_dims_stage *dims);
//
ocp_qp_in_stage *ocp_qp_in_stage_assign(void *config, ocp_qp_dims_stage *dims, void *raw_memory);
//
int ocp_qp_out_calculate_size(void *config, ocp_qp_dims *dims);
//
ocp_qp_out *ocp_qp_out_assign(void *config, ocp_qp_dims *dims, void *raw_memory);
//
int ocp_qp_res_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_res *ocp_qp_res_assign(ocp_qp_dims *dims, void *raw_memory);
//
int ocp_qp_res_workspace_calculate_size(ocp_qp_dims *dims);
//
ocp_qp_res_ws *ocp_qp_res_workspace_assign(ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_res_compute(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_qp_res *qp_res, ocp_qp_res_ws *res_ws);
//
void ocp_qp_res_compute_nrm_inf(ocp_qp_res *qp_res, double res[4]);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_QP_OCP_QP_COMMON_H_

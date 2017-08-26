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

#ifndef ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_
#define ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_

#ifdef __cplusplus
extern "C" {
#endif



// XXX comment this in for HPIPM condensing, out for ACADOS condensing
#define HPIPM_COND



#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

// giaf condensing
#if defined(HPIPM_COND)



#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"



// define maximum number of variables and constraints in qpoases XXX does not work !
//#define QPOASES_NVMAX 200
//#define QPOASES_NCMAX 600
//#define __EXTERNAL_DIMENSIONS__



// struct of arguments to the solver
typedef struct ocp_qp_condensing_qpoases_args_ {
	double cputime; // maximum cpu time in seconds
	int nwsr; // maximum number of working set recalculations
	int warm_start; // warm start with dual_sol in memory
} ocp_qp_condensing_qpoases_args;



// struct of the solver memory
typedef struct ocp_qp_condensing_qpoases_memory_ {
	struct d_ocp_qp *qp;
	struct d_ocp_qp_sol *qp_sol;
	struct d_dense_qp *qpd;
	struct d_dense_qp_sol *qpd_sol;
	struct d_cond_qp_ocp2dense_workspace *cond_workspace;
	struct d_strmat *sR;
	double **hlam_lb;
	double **hlam_ub;
	double **hlam_lg;
	double **hlam_ug;
	double *H;
	double *R;
	double *g;
	double *A;
	double *b;
	double *d_lb0;
	double *d_ub0;
	double *d_lb;
	double *d_ub;
	double *C;
	double *d_lg;
	double *d_ug;
	int *idxb;
	double *prim_sol;
	double *dual_sol;
	void *QPB; // XXX cast to QProblemB to use !!!
	void *QP; // XXX cast to QProblem to use !!!
	double inf_norm_res[5];
	double cputime; // requried cpu time
	int nwsr; // performed number of working set recalculations
} ocp_qp_condensing_qpoases_memory;



//
int ocp_qp_condensing_qpoases_calculate_workspace_size(ocp_qp_in *qp_in, ocp_qp_condensing_qpoases_args *args);
//
int ocp_qp_condensing_qpoases_calculate_memory_size(ocp_qp_in *qp_in, ocp_qp_condensing_qpoases_args *args);
//
void ocp_qp_condensing_qpoases_create_memory(ocp_qp_in *qp_in, ocp_qp_condensing_qpoases_args *args, ocp_qp_condensing_qpoases_memory *qpoases_memory, void *memory);
//
int_t ocp_qp_condensing_qpoases(ocp_qp_in *input, ocp_qp_out *output,
    void *args, void *mem, void *work);


// XXX remove !!!!!
void ocp_qp_condensing_qpoases_initialize(ocp_qp_in *qp_in, void *args_, void *mem_, void **work);
void ocp_qp_condensing_qpoases_destroy(void *mem, void *work);



// robin condensing
#else



typedef struct {
    real_t dummy;
} ocp_qp_condensing_qpoases_args;

int_t ocp_qp_condensing_qpoases(ocp_qp_in *input, ocp_qp_out *output,
    void *args, void *mem, void *work);

int_t ocp_qp_condensing_qpoases_workspace_size(ocp_qp_in *input,
    ocp_qp_condensing_qpoases_args *args);

void ocp_qp_condensing_qpoases_initialize(ocp_qp_in *qp_in, void *args_, void *mem_, void **work);
void ocp_qp_condensing_qpoases_destroy(void *mem, void *work);



// giaf vs robin condensing
#endif



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif  // ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_

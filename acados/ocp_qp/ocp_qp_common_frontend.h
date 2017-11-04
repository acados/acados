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

#ifndef ACADOS_OCP_QP_OCP_QP_COMMON_FRONTEND_H_
#define ACADOS_OCP_QP_OCP_QP_COMMON_FRONTEND_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"


typedef struct {
    int N;
    int *nx;
    int *nu;
    int *nb;
    int *nc;
    double **A;
    double **B;
    double **b;
    double **Q;
    double **S;
    double **R;
    double **q;
    double **r;
    int **idxb;
    double **lb;
    double **ub;
    double **Cx;
    double **Cu;
    double **lc;
    double **uc;
} col_maj_ocp_qp_in;


typedef struct {
    double **x;
    double **u;
    double **pi;
    double **lam;
} col_maj_ocp_qp_out;


//
int col_maj_ocp_qp_in_calculate_size(ocp_qp_dims *dims);
//
char *assign_col_maj_ocp_qp_in(ocp_qp_dims *dims, col_maj_ocp_qp_in **qp_in, void *ptr);
//
int col_maj_ocp_qp_out_calculate_size(ocp_qp_dims *dims);
//
char *assign_col_maj_ocp_qp_out(ocp_qp_dims *dims, col_maj_ocp_qp_out **qp_out, void *ptr);
//
void convert_col_maj_ocp_qp_out(ocp_qp_dims *dims, ocp_qp_out *qp_out, col_maj_ocp_qp_out *cm_qp_out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_QP_OCP_QP_COMMON_FRONTEND_H_

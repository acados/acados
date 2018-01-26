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
} colmaj_ocp_qp_in;


typedef struct {
    double **x;
    double **u;
    double **pi;
    double **lam;
} colmaj_ocp_qp_out;


typedef struct {
    double **res_r;
    double **res_q;
    double **res_ls;
    double **res_us;
    double **res_b;
    double **res_d_lb;
    double **res_d_ub;
    double **res_d_lg;
    double **res_d_ug;
    double **res_d_ls;
    double **res_d_us;
    double **res_m_lb;
    double **res_m_ub;
    double **res_m_lg;
    double **res_m_ug;
    double **res_m_ls;
    double **res_m_us;
    double res_nrm_inf[4];
} colmaj_ocp_qp_res;


//
int colmaj_ocp_qp_in_calculate_size(ocp_qp_dims *dims);
//
char *assign_colmaj_ocp_qp_in(ocp_qp_dims *dims, colmaj_ocp_qp_in **qp_in, void *ptr);
//
int colmaj_ocp_qp_out_calculate_size(ocp_qp_dims *dims);
//
char *assign_colmaj_ocp_qp_out(ocp_qp_dims *dims, colmaj_ocp_qp_out **qp_out, void *ptr);
//
int colmaj_ocp_qp_res_calculate_size(ocp_qp_dims *dims);
//
char *assign_colmaj_ocp_qp_res(ocp_qp_dims *dims, colmaj_ocp_qp_res **qp_res, void *ptr);
//
void convert_colmaj_to_ocp_qp_in(colmaj_ocp_qp_in *cm_qp_in, ocp_qp_in *qp_in);
//
void convert_ocp_qp_out_to_colmaj(ocp_qp_out *qp_out, colmaj_ocp_qp_out *cm_qp_out);
//
void convert_ocp_qp_res_to_colmaj(ocp_qp_res *qp_res, colmaj_ocp_qp_res *cm_qp_res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_QP_OCP_QP_COMMON_FRONTEND_H_

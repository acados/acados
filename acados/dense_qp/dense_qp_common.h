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


#ifndef ACADOS_DENSE_QP_DENSE_QP_COMMON_H_
#define ACADOS_DENSE_QP_DENSE_QP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// hpipm
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
// acados
#include "acados/utils/types.h"


typedef struct d_dense_qp dense_qp_in;
typedef struct d_dense_qp_sol dense_qp_out;


//
int dense_qp_in_calculate_size(int nv, int ne, int nb, int ng, int ns);
//
char *assign_dense_qp_in(int nv, int ne, int nb, int ng, int ns, dense_qp_in **qp_in, void *ptr);
//
int dense_qp_out_calculate_size(int nv, int ne, int nb, int ng, int ns);
//
char *assign_dense_qp_out(int nv, int ne, int nb, int ng, int ns, dense_qp_out **qp_out, void *ptr);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_DENSE_QP_DENSE_QP_COMMON_H_

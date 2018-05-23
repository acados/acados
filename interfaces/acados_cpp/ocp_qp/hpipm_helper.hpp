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

#ifndef INTERFACES_ACADOS_CPP_OCP_QP_HPIPM_HELPER_HPP_
#define INTERFACES_ACADOS_CPP_OCP_QP_HPIPM_HELPER_HPP_

#define num_rows_Q(stage, dim) (dim->nx[stage])

#define num_cols_Q(stage, dim) (dim->nx[stage])

#define num_rows_S(stage, dim) (dim->nu[stage])

#define num_cols_S(stage, dim) (dim->nx[stage])

#define num_rows_R(stage, dim) (dim->nu[stage])

#define num_cols_R(stage, dim) (dim->nu[stage])

#define num_elems_q(stage, dim) (dim->nx[stage])

#define num_elems_r(stage, dim) (dim->nu[stage])

#define num_rows_A(stage, dim) (dim->nx[stage + 1])

#define num_cols_A(stage, dim) (dim->nx[stage])

#define num_rows_B(stage, dim) (dim->nx[stage + 1])

#define num_cols_B(stage, dim) (dim->nu[stage])

#define num_elems_b(stage, dim) (dim->nx[stage + 1])

#define num_elems_lbx(stage, dim) (dim->nbx[stage])

#define num_elems_lbu(stage, dim) (dim->nbu[stage])

#define num_elems_ubx(stage, dim) (dim->nbx[stage])

#define num_elems_ubu(stage, dim) (dim->nbu[stage])

#define num_rows_C(stage, dim) (dim->ng[stage])

#define num_cols_C(stage, dim) (dim->nx[stage])

#define num_rows_D(stage, dim) (dim->ng[stage])

#define num_cols_D(stage, dim) (dim->nu[stage])

#define num_elems_lg(stage, dim) (dim->ng[stage])

#define num_elems_ug(stage, dim) (dim->ng[stage])

#endif  // INTERFACES_ACADOS_CPP_OCP_QP_HPIPM_HELPER_HPP_

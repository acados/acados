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

#ifndef ACADOS_SIM_CASADI_WRAPPER_H_
#define ACADOS_SIM_CASADI_WRAPPER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

void vde_fun(const int_t nx, const int_t nu, const real_t *in, real_t *out, casadi_function_t vde);

void jac_fun(const int_t nx, const real_t *in, real_t *out, casadi_function_t jac);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_CASADI_WRAPPER_H_

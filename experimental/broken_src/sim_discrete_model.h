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

#ifndef ACADOS_SIM_SIM_DISCRETE_MODEL_H_
#define ACADOS_SIM_SIM_DISCRETE_MODEL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_common.h"

int sim_discrete_model(const sim_in *, sim_out *, void *, void *, void *);

//void vde_adj_fun(const int nx, const int nu, const double *in, double *out,
//                 casadi_function_t vde_hess);

//void vde_hess_fun(const int nx, const int nu, const double *in, double *out,
//                  casadi_function_t vde_hess);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_DISCRETE_MODEL_H_

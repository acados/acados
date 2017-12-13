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

#ifndef ACADOS_SIM_SIM_COLLOCATION_H_
#define ACADOS_SIM_SIM_COLLOCATION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

enum Newton_type_collocation { exact = 0, simplified_in, simplified_inis };

typedef struct {
    enum Newton_type_collocation type;
    real_t *eig;
    real_t *low_tria;
    bool single;
    bool freeze;

    real_t *transf1;
    real_t *transf2;

    real_t *transf1_T;
    real_t *transf2_T;
} Newton_scheme;

real_t LU_system_solve(real_t *const A, real_t *const b, int *const perm,
                       int dim, int dim2);

void get_Gauss_nodes(const int_t num_stages, real_t *nodes);

void read_Gauss_simplified(const int_t num_stages, Newton_scheme *scheme);

void create_Butcher_table(const int_t num_stages, const real_t *nodes,
                          real_t *b, real_t *A);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_COLLOCATION_H_

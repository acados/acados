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

typedef struct simplified_form_ {
    real_t *eig;

    real_t *transf1;
    real_t *transf2;

    real_t *transf1_T;
    real_t *transf2_T;
} simplified_form;

typedef struct single_form {
    real_t eig;
    real_t *low_tria;

    real_t *transf1;
    real_t *transf2;

    real_t *transf1_T;
    real_t *transf2_T;
} single_form;

void get_Gauss_nodes(const int_t num_stages, real_t *nodes);

void create_Butcher_table(const int_t num_stages, const real_t *nodes,
        real_t *b, real_t *A);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_COLLOCATION_H_

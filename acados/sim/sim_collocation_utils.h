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

#ifndef ACADOS_SIM_SIM_COLLOCATION_UTILS_H_
#define ACADOS_SIM_SIM_COLLOCATION_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"



enum Newton_type_collocation
{
    exact = 0,
    simplified_in,
    simplified_inis
};



typedef struct
{
    enum Newton_type_collocation type;
    double *eig;
    double *low_tria;
    bool single;
    bool freeze;

    double *transf1;
    double *transf2;

    double *transf1_T;
    double *transf2_T;
} Newton_scheme;



//
int gauss_nodes_work_calculate_size(int ns);
//
void gauss_nodes(int ns, double *nodes, void *raw_memory);
//
int gauss_simplified_work_calculate_size(int ns);
//
void gauss_simplified(int ns, Newton_scheme *scheme, void *work);
//
int butcher_table_work_calculate_size(int ns);
//
void butcher_table(int ns, double *nodes, double *b, double *A, void *work);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_COLLOCATION_UTILS_H_

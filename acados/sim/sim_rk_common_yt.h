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

#ifndef ACADOS_SIM_SIM_RK_COMMON_YT_H_
#define ACADOS_SIM_SIM_RK_COMMON_YT_H_

#include "acados/utils/types.h"

typedef struct {
    int num_stages;
    int newton_iter;
    double *A_mat;
    double *c_vec;
    double *b_vec;  
} sim_RK_opts;

int_t sim_RK_opts_calculate_size(int_t ns);

char *assign_sim_RK_opts(int_t ns, sim_RK_opts **opts, void *ptr);

sim_RK_opts *create_sim_RK_opts(int_t ns);

#endif  // ACADOS_SIM_SIM_RK_COMMON_YT_H_

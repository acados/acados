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

// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_gnsf.h"
//#include "acados/sim/sim_erk_integrator.h"
//#include "acados/sim/sim_lifted_irk_integrator.h"



void print_gnsf_dims( gnsf_dims dims)
{
    printf("\n");
    printf("nx  = %d \n", dims.nx);
    printf("nu  = %d \n", dims.nu);
    printf("nz  = %d \n", dims.nz);
    printf("nx1 = %d \n", dims.nx1);
    printf("nx2 = %d \n", dims.nx2);
    printf("n_out = %d \n", dims.n_out);
    printf("n_steps  = %d \n", dims.num_steps);
    printf("n_stages = %d \n", dims.num_stages);

    printf("\n");
}
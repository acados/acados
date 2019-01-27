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

#include "acados/sim/sim_discrete_model.h"

#include <stdlib.h>

#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/sim/sim_common.h"

int sim_discrete_model(const sim_in *in, sim_out *out, void *args, void *mem, void *work)
{
    int nx = in->nx;
    int nu = in->nu;
    double *in_array, *out_array;
	// TODO use memory instead of calloc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    in_array = calloc(nx + nu, sizeof(double));
    out_array = calloc(nx + nx * nu, sizeof(double));
    for (int i = 0; i < nx; i++) in_array[i] = in->x[i];
    for (int i = 0; i < nu; i++) in_array[nx + i] = in->u[i];
    discrete_model_fun(in->nx, in->nu, in_array, out_array, in->discrete_model);
    for (int i = 0; i < nx; i++) out->xn[i] = out_array[i];
    for (int i = 0; i < nx * (nx + nu); i++) out->S_forw[i] = out_array[nx + i];
	free(in_array);
	free(out_array);
    return 0;
}

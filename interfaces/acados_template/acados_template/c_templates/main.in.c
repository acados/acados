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
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_{{ra.model_name}}.h"

int main() {

    int status = 0;
    status = acados_create();

    if (status) { 
        printf("acados_create() returned status %d. Exiting.\n", status); 
        exit(1); }

    status = acados_solve();

    if (status) { 
        printf("acados_solve() returned status %d. Exiting.\n", status); 
        exit(1); }

    status = acados_free();

    if (status) { 
        printf("acados_free() returned status %d. Exiting.\n", status); 
        exit(1); }

    return status;
}

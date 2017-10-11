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

#include "acados/ocp_nlp/allocate_ocp_nlp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

void allocate_ocp_nlp_in(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *ng, ocp_nlp_in *const nlp_in) {
    // TODO(nielsvd): implement.
}

void free_ocp_nlp_in(ocp_nlp_in *const nlp) {
    // TODO(nielsvd): implement.
}

void allocate_ocp_nlp_out(ocp_nlp_in *const in, ocp_nlp_out *out) {
    // TODO(nielsvd): implement.
}

void free_ocp_nlp_out(int_t N, ocp_nlp_out *out) {
    // TODO(nielsvd): implement.
}

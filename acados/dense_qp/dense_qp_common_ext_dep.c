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

// external
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
// acados
#include "acados/dense_qp/dense_qp_common_ext_dep.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"



dense_qp_in *create_dense_qp_in(int nv, int ne, int nb, int ng, int ns)
{
    int size = dense_qp_in_calculate_size(nv, ne, nb, ng, ns);
    void *ptr = malloc(size);
    dense_qp_in *qp_in = assign_dense_qp_in(nv, ne, nb, ng, ns, ptr);
    return qp_in;
}



dense_qp_out *create_dense_qp_out(int nv, int ne, int nb, int ng, int ns)
{
    int size = dense_qp_out_calculate_size(nv, ne, nb, ng, ns);
    void *ptr = malloc(size);
    dense_qp_out *qp_out = assign_dense_qp_out(nv, ne, nb, ng, ns, ptr);
    return qp_out;
}



void print_dense_qp_in(dense_qp_in *qp_in)
{
    int nv = qp_in->nv;

    printf("H =\n");
    d_print_strmat(nv, nv, qp_in->Hg, 0, 0);
    // TODO(dimitris): print all data
}
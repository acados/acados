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
// // blasfeo
// #include "blasfeo_target.h"
// #include "blasfeo_common.h"
// #include "blasfeo_d_aux.h"
// #include "blasfeo_d_blas.h"
// #include "blasfeo_v_aux_ext_dep.h"
// // hpipm
// #include "hpipm_d_ocp_qp.h"
// #include "hpipm_d_ocp_qp_sol.h"
// #include "hpipm_d_ocp_qp_ipm.h"
// acados
#include "acados/ocp_qp/ocp_qp_condensing_hpipm_ext_dep.h"
#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_common.h"
// #include "utils/types.h"


ocp_qp_condensing_hpipm_args *ocp_qp_condensing_hpipm_create_arguments(ocp_qp_dims *dims) {

    int memory_size = ocp_qp_condensing_hpipm_calculate_args_size(dims);
    void *ptr = malloc(memory_size);
    ocp_qp_condensing_hpipm_args *args = ocp_qp_condensing_hpipm_assign_args(dims, ptr);
    ocp_qp_condensing_hpipm_initialize_default_args(args);
    return args;
}



ocp_qp_condensing_hpipm_memory *ocp_qp_condensing_hpipm_create_memory(ocp_qp_dims *dims, void *args_) {
    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;

    ocp_qp_condensing_hpipm_memory *mem;

    int memory_size = ocp_qp_condensing_hpipm_calculate_memory_size(dims, args);
    void *ptr = malloc(memory_size);
    char *ptr_end = ocp_qp_condensing_hpipm_assign_memory(dims, args, (void **) &mem, ptr);
    assert((char *) ptr + memory_size >= ptr_end); (void) ptr_end;

    return mem;
}



void ocp_qp_condensing_hpipm_initialize(ocp_qp_in *qp_in, void *args_, void **mem, void **work) {
    ocp_qp_condensing_hpipm_args *args = (ocp_qp_condensing_hpipm_args *) args_;

    *mem = ocp_qp_condensing_hpipm_create_memory(qp_in->size, args);

}



void ocp_qp_condensing_hpipm_destroy(void *mem, void *work) {
    free(mem);
    free(work);
}

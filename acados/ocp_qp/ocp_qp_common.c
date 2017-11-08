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

// hpipm
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp_size.h"
// acados
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"

//#include "acados/ocp_qp/ocp_qp_hpipm.h"
//#ifdef ACADOS_WITH_HPMPC
//#include "acados/ocp_qp/ocp_qp_hpmpc.h"
//#endif
//#ifdef OOQP
//#include "acados/ocp_qp/ocp_qp_ooqp.h"
//#endif
//#include "acados/ocp_qp/ocp_qp_qpdunes.h"
//#include "acados/utils/types.h"



int ocp_qp_in_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_in);
    size += d_memsize_ocp_qp(dims);
    size += sizeof(ocp_qp_dims);
    size += d_memsize_ocp_qp_size(dims->N);
    return size;
}



char *assign_ocp_qp_in(ocp_qp_dims *dims, ocp_qp_in **qp_in, void *mem)
{
    // char pointer
    char *c_ptr = (char *) mem;

    *qp_in = (ocp_qp_in *) c_ptr;
    c_ptr += sizeof(ocp_qp_in);

    ocp_qp_dims *dims_copy = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    d_create_ocp_qp(dims, *qp_in, c_ptr);
    c_ptr += d_memsize_ocp_qp(dims);

    d_create_ocp_qp_size(dims->N, dims_copy, c_ptr);
    c_ptr += d_memsize_ocp_qp_size(dims->N);

    dims_copy->N = dims->N;

    for (int ii = 0; ii < dims->N+1; ii++)
        {
        dims_copy->nx[ii] = dims->nx[ii];
        dims_copy->nu[ii] = dims->nu[ii];
        dims_copy->nb[ii] = dims->nb[ii];
        dims_copy->ng[ii] = dims->ng[ii];
        dims_copy->ns[ii] = dims->ns[ii];
        dims_copy->nbu[ii] = dims->nbu[ii];
        dims_copy->nbx[ii] = dims->nbx[ii];
        }

    (*qp_in)->size = dims_copy;

    return c_ptr;
}



int ocp_qp_out_calculate_size(ocp_qp_dims *dims)
{
    int size = sizeof(ocp_qp_out);
    size += d_memsize_ocp_qp_sol(dims);

    return size;
}



char *assign_ocp_qp_out(ocp_qp_dims *dims, ocp_qp_out **qp_out, void *mem)
{
    // char pointer
    char *c_ptr = (char *) mem;

    *qp_out = (ocp_qp_out *) c_ptr;
    c_ptr += sizeof(ocp_qp_out);

    d_create_ocp_qp_sol(dims, *qp_out, c_ptr);
    c_ptr += d_memsize_ocp_qp_sol(dims);

    return c_ptr;
}




void form_nbu_nbx_rev(int N, int *nbu, int *nbx, int *nb, int* nx, int *nu, int **idxb_rev)
{
    for (int ii = 0; ii < N+1; ii++)
    {
        nbu[ii] = 0;
        nbx[ii] = 0;
        for (int jj = 0; jj < nb[ii]; jj++)
        {
            if (idxb_rev[ii][jj] < nx[ii])
            {
                nbx[ii]++;
            }
            else
            {
                nbu[ii]++;
            }
        }
    }
}


new_ocp_qp_solver initialize_ocp_qp_solver(qp_solver_t qp_solver_name)
{
    new_ocp_qp_solver qp_solver;

    qp_solver.qp_in = NULL;
    qp_solver.qp_out = NULL;
    qp_solver.args = NULL;
    qp_solver.mem = NULL;
    qp_solver.work = NULL;

    switch (qp_solver_name)
    {
        case HPIPM:
            qp_solver.calculate_args_size = &ocp_qp_hpipm_calculate_args_size;
            qp_solver.assign_args = &ocp_qp_hpipm_assign_args;
            // ...
            break;
        case CONDENSING_HPIPM:
            qp_solver.calculate_args_size = &ocp_qp_condensing_hpipm_calculate_args_size;
            qp_solver.assign_args = &ocp_qp_condensing_hpipm_assign_args;
            break;
        case CONDENSING_QPOASES:
            qp_solver.calculate_args_size = &ocp_qp_condensing_qpoases_calculate_args_size;
            qp_solver.assign_args = &ocp_qp_condensing_qpoases_assign_args;
            qp_solver.initialize_default_args = &ocp_qp_condensing_qpoases_initialize_default_args;
            break;
        default:
            qp_solver.calculate_args_size = NULL;
    }
    return qp_solver;
}

// NOTE(dimitris): maybe better do switch only in one function, as above
// int ocp_qp_solver_calculate_args_size(ocp_qp_dims *dims, qp_solver_t qp_solver)
// {
//     int size = 0;

//     switch (qp_solver)
//     {
//         case HPIPM:
//             size = ocp_qp_hpipm_calculate_args_size(dims);
//         case CONDENSING_HPIPM:
//             size = ocp_qp_condensing_hpipm_calculate_args_size(dims);
//         case CONDENSING_QPOASES:
//             size = ocp_qp_condensing_qpoases_calculate_args_size(dims);
//         default:
//             size = -1;
//     }
//     return size;
// }

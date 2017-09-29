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

#ifndef ACADOS_OCP_NLP_OCP_NLP_SQP_H_
#define ACADOS_OCP_NLP_OCP_NLP_SQP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/utils/types.h"

int_t ocp_nlp_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out,
                     void *nlp_args, void *nlp_mem, void *work);

void ocp_nlp_sqp_create_memory(const ocp_nlp_in *in, void *args_,
                                  void *memory_);

void ocp_nlp_sqp_free_memory(void *memory_);

int_t ocp_nlp_sqp_calculate_workspace_size(const ocp_nlp_in *in,
                                              void *args_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_SQP_H_

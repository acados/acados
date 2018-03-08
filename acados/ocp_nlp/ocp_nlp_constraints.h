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

#ifndef ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_H_
#define ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
// acados
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"
#include "acados/ocp_qp/ocp_qp_common.h"



/************************************************
* dims
************************************************/

typedef struct
{
    int nx;
    int nu;
    int nb;  // nbx + nbu
    int nbx;
    int nbu;
    int ng;  // number of general linear constraints
    int ns;  // number of soft constraints
} ocp_nlp_constraints_dims;

//
int ocp_nlp_constraints_dims_calculate_size();
//
ocp_nlp_constraints_dims *ocp_nlp_constraints_dims_assign(void *raw_memory);



/************************************************
* constraints
************************************************/

typedef struct
{
	ocp_nlp_constraints_dims *dims;
    int *idxb;
	struct blasfeo_dvec d;
	struct blasfeo_dmat DCt;
}
ocp_nlp_constraints_model;

//
int ocp_nlp_constraints_model_calculate_size(void *config, ocp_nlp_constraints_dims *dims);
//
void *ocp_nlp_constraints_model_assign(void *config, ocp_nlp_constraints_dims *dims, void *raw_memory);
//
void ocp_nlp_constraints_initialize_qp(void *config, ocp_nlp_constraints_dims *dims, ocp_nlp_constraints_model *model, ocp_qp_in_stage *qp_in_stage, void *mem, void *work); // TODO mem and work if needed



/************************************************
* config
************************************************/

typedef struct
{
	int (*model_calculate_size) (void *config, ocp_nlp_constraints_dims *dims);
	void *(*model_assign) (void *config, ocp_nlp_constraints_dims *dims, void *raw_memory);
	void (*initialize_qp) (void *config, ocp_nlp_constraints_dims *dims, ocp_nlp_constraints_model *model, ocp_qp_in_stage *qp_in_stage, void *mem, void *work);
	void (*config_initialize_default) (void *config);
} ocp_nlp_constraints_config;

//
int ocp_nlp_constraints_config_calculate_size();
//
ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory);
//
void ocp_nlp_constraints_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_H_

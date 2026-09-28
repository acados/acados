/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


/// \addtogroup ocp_nlp
/// @{
/// \addtogroup ocp_nlp_reg
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_REG_NOREG_H_
#define ACADOS_OCP_NLP_OCP_NLP_REG_NOREG_H_

#ifdef __cplusplus
extern "C" {
#endif



// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"



/************************************************
 * dims
 ************************************************/

// use the functions in ocp_nlp_reg_common


/************************************************
 * options
 ************************************************/

//
acados_size_t ocp_nlp_reg_noreg_opts_calculate_size(void);
//
void *ocp_nlp_reg_noreg_opts_assign(void *raw_memory);
//
void ocp_nlp_reg_noreg_opts_initialize_default(void *config_, ocp_nlp_reg_dims *dims, void *opts_);
//
void ocp_nlp_reg_noreg_opts_set(void *config_, void *opts_, const char *field, void* value);



/************************************************
 * memory
 ************************************************/
//
acados_size_t ocp_nlp_reg_noreg_memory_calculate_size(void *config, ocp_nlp_reg_dims *dims, void *opts);
//
void *ocp_nlp_reg_noreg_memory_assign(void *config, ocp_nlp_reg_dims *dims, void *opts, void *raw_memory);

/************************************************
 * workspace
 ************************************************/

// not needed

/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_reg_noreg_config_initialize_default(ocp_nlp_reg_config *config);



#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_REG_NOREG_H_

/// @}
/// @}

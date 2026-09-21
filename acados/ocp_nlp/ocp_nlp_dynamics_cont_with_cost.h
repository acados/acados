/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 */

#ifndef ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_CONT_WITH_COST_H_
#define ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_CONT_WITH_COST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"

void ocp_nlp_dynamics_cont_with_cost_config_initialize_default(void *config, int stage);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_CONT_WITH_COST_H_
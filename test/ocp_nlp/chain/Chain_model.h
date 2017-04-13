#ifndef TEST_OCP_NLP_CHAIN_CHAIN_MODEL_H_
#define TEST_OCP_NLP_CHAIN_CHAIN_MODEL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

void VDE_fun_nm2(const real_t* in, real_t* out,
    int (*vde)(const real_t**, real_t**, int*, real_t*, int));
void VDE_fun_nm3(const real_t* in, real_t* out,
    int (*vde)(const real_t**, real_t**, int*, real_t*, int));
void VDE_fun_nm4(const real_t* in, real_t* out,
    int (*vde)(const real_t**, real_t**, int*, real_t*, int));

void jac_fun_nm2(const real_t* in, real_t* out);
void jac_fun_nm3(const real_t* in, real_t* out);
void jac_fun_nm4(const real_t* in, real_t* out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // TEST_OCP_NLP_CHAIN_CHAIN_MODEL_H_

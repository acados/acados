#ifndef TEST_OCP_NLP_CHAIN_CHAIN_MODEL_H_
#define TEST_OCP_NLP_CHAIN_CHAIN_MODEL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

int vde_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int vde_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int vde_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

int jac_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int jac_chain_nm3(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int jac_chain_nm4(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // TEST_OCP_NLP_CHAIN_CHAIN_MODEL_H_

#ifndef EXAMPLES_CASADI_CHAIN_CHAIN_MODEL_H_
#define EXAMPLES_CASADI_CHAIN_CHAIN_MODEL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

void VDE_fun_nm2(const real_t* in, real_t* out);
void VDE_fun_nm3(const real_t* in, real_t* out);
void VDE_fun_nm4(const real_t* in, real_t* out);
void VDE_fun_nm5(const real_t* in, real_t* out);
//void VDE_fun_nm6(const real_t* in, real_t* out);
//void VDE_fun_nm7(const real_t* in, real_t* out);
//void VDE_fun_nm8(const real_t* in, real_t* out);
//void VDE_fun_nm9(const real_t* in, real_t* out);

void jac_fun_nm2(const real_t* in, real_t* out);
void jac_fun_nm3(const real_t* in, real_t* out);
void jac_fun_nm4(const real_t* in, real_t* out);
void jac_fun_nm5(const real_t* in, real_t* out);
//void jac_fun_nm6(const real_t* in, real_t* out);
//void jac_fun_nm7(const real_t* in, real_t* out);
//void jac_fun_nm8(const real_t* in, real_t* out);
//void jac_fun_nm9(const real_t* in, real_t* out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_CASADI_CHAIN_CHAIN_MODEL_H_

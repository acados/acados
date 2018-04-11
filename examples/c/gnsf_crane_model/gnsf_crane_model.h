#ifndef EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_
#define EXAMPLES_C_GNSF_CRANE_MODEL_GNSF_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

// used to import integers & double matrices
int get_ints_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int get_matrices_fun(const double** arg, double** res, int* iw, double* w, void *mem);
// int But_KK_YY_ZZ_LO_fun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

// phi_fun_jac_y
int        phi_fun_jac_y(const double** arg, double** res, int* iw, double* w, void *mem);
int        phi_fun_jac_y_work(int *, int *, int *, int *);
const int *phi_fun_jac_y_sparsity_in(int);
const int *phi_fun_jac_y_sparsity_out(int);
int        phi_fun_jac_y_n_in();
int        phi_fun_jac_y_n_out();

// phi_jac_y_uhat
int        phi_jac_y_uhat(const double** arg, double** res, int* iw, double* w, void *mem);
int        phi_jac_y_uhat_work(int *, int *, int *, int *);
const int *phi_jac_y_uhat_sparsity_in(int);
const int *phi_jac_y_uhat_sparsity_out(int);
int        phi_jac_y_uhat_n_in();
int        phi_jac_y_uhat_n_out();

// f_lo_fun_jac_x1k1uz
int        f_lo_fun_jac_x1k1uz(const double** arg, double** res, int* iw, double* w, void *mem);
int        f_lo_fun_jac_x1k1uz_work(int *, int *, int *, int *);
const int *f_lo_fun_jac_x1k1uz_sparsity_in(int);
const int *f_lo_fun_jac_x1k1uz_sparsity_out(int);
int        f_lo_fun_jac_x1k1uz_n_in();
int        f_lo_fun_jac_x1k1uz_n_out();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // EXAMPLES_C_CRANE_MODEL_CRANE_MODEL_H_
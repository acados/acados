#include <math.h>
#include <casadi/casadi.hpp>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

bool  libexternal_ode_casadi_has_derivative(void);
const char* libexternal_ode_casadi_name(void);

int libexternal_ode_casadi(const casadi_real** arg,
                                                   double** res,
                                                   casadi_int* iw,
                                                   casadi_real* w,
                                                   void* mem);

// IN
casadi_int libexternal_ode_casadi_n_in(void);
const char* libexternal_ode_casadi_name_in(casadi_int i);
const casadi_int* libexternal_ode_casadi_sparsity_in(casadi_int i);

// OUT
casadi_int libexternal_ode_casadi_n_out(void);
const char* libexternal_ode_casadi_name_out(casadi_int i);
const casadi_int* libexternal_ode_casadi_sparsity_out(casadi_int i);

int libexternal_ode_casadi_work(casadi_int *sz_arg,
                                               casadi_int* sz_res,
                                               casadi_int *sz_iw,
                                               casadi_int *sz_w);

#ifdef __cplusplus
} /* extern "C" */
#endif


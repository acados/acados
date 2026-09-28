/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#ifndef {{ name }}_p_global_PRECOMPUTE_FUN
#define {{ name }}_p_global_PRECOMPUTE_FUN

#ifdef __cplusplus
extern "C" {
#endif

int {{ name }}_p_global_precompute_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ name }}_p_global_precompute_fun_work(int *, int *, int *, int *);
const int *{{ name }}_p_global_precompute_fun_sparsity_in(int);
const int *{{ name }}_p_global_precompute_fun_sparsity_out(int);
int {{ name }}_p_global_precompute_fun_n_in(void);
int {{ name }}_p_global_precompute_fun_n_out(void);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ name }}_p_global_PRECOMPUTE_FUN

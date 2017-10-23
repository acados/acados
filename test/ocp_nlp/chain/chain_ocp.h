#ifndef TEST_OCP_NLP_CHAIN_CHAIN_OCP_H_
#define TEST_OCP_NLP_CHAIN_CHAIN_OCP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

int ls_cost_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int ls_cost_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int ls_cost_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);

int ls_cost_nm2_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm3_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int ls_cost_nm4_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);

const int *ls_cost_nm2_sparsity_out(int i);
const int *ls_cost_nm3_sparsity_out(int i);
const int *ls_cost_nm4_sparsity_out(int i);

int ls_costN_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int ls_costN_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int ls_costN_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);

int ls_costN_nm2_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm3_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int ls_costN_nm4_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);

const int *ls_costN_nm2_sparsity_out(int i);
const int *ls_costN_nm3_sparsity_out(int i);
const int *ls_costN_nm4_sparsity_out(int i);

int pathcon_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int pathcon_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int pathcon_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);

int pathcon_nm2_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int pathcon_nm3_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int pathcon_nm4_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);

const int *pathcon_nm2_sparsity_out(int i);
const int *pathcon_nm3_sparsity_out(int i);
const int *pathcon_nm4_sparsity_out(int i);

int pathconN_nm2(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int pathconN_nm3(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);
int pathconN_nm4(const real_t **arg, real_t **res, int *iw, real_t *w, int mem);

int pathconN_nm2_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int pathconN_nm3_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);
int pathconN_nm4_work(int *sz_arg, int *sz_res, int *sz_iw, int *sz_w);

const int *pathconN_nm2_sparsity_out(int i);
const int *pathconN_nm3_sparsity_out(int i);
const int *pathconN_nm4_sparsity_out(int i);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // TEST_OCP_NLP_CHAIN_CHAIN_OCP_H_

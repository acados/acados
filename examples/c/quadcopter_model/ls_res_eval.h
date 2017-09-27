#ifndef LS_RES_EVAL_H_
#define LS_RES_EVAL_H_

#include "acados/utils/types.h"
// #include "acados_kitty_interface.h"

void ls_res_eval(const real_t* in, real_t* out,
                      int (*ls_res)(const real_t**, real_t**, int*, real_t*,
                                 int)) {
#endif  // LS_RES_EVAL_H_

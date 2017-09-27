#ifndef LS_RES_EVAL_H_
#define LS_RES_EVAL_H_

#include "acados/utils/types.h"

void ls_res_eval_quadcopter(const real_t* in, real_t* out,
                      int (*ls_res)(const real_t**, real_t**, int*, real_t*, int));

void ls_res_eval_end_quadcopter(const real_t* in, real_t* out,
                      int (*ls_res)(const real_t**, real_t**, int*, real_t*, int));

#endif  // LS_RES_EVAL_H_

#ifndef RESIDUAL_X_EVAL_WRAPPER_H_
#define RESIDUAL_X_EVAL_WRAPPER_H_

#include "acados/utils/types.h"
#include "acados_kitty_interface.h"

void residual_x_eval_wrapper(nmpc_data* nmpc_data, const real_t* in, real_t* out,  void *internal_mem);

#endif  // RESIDUAL_X_EVAL_WRAPPER_H_

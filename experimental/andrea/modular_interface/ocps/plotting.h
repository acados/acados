#ifndef KITTY_ACADOS_PLOTTING_H_
#define KITTY_ACADOS_PLOTTING_H_
#include <stdio.h>
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/utils/tools.h"
void plot_states_controls(real_t *w, real_t T, int_t NN, int_t NX, int_t NU, FILE *gnuplotPipe);
#endif

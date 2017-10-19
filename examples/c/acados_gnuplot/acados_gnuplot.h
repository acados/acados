#ifndef EXAMPLES_C_ACADOS_GNUPLOT_ACADOS_GNUPLOT_H_
#define EXAMPLES_C_ACADOS_GNUPLOT_ACADOS_GNUPLOT_H_

#include "acados/utils/types.h"

void acados_gnuplot(real_t **data, int_t n_data, real_t T,
    int_t N, char **labels, int_t layout_x, int_t layout_y);

#endif  // EXAMPLES_C_ACADOS_GNUPLOT_ACADOS_GNUPLOT_H_

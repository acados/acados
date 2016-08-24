#ifndef ACADOS_PRINT_H_
#define ACADOS_PRINT_H_

#include "acados/types.h"

void print_matrix(char *file_name, const real_t *matrix, const int_t nrows,
    const int_t ncols);

void print_array(char *file_name, real_t *array, int_t size);

#endif  // ACADOS_PRINT_H_

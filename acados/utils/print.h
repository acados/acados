#ifndef ACADOS_UTILS_PRINT_H_
#define ACADOS_UTILS_PRINT_H_

#include "acados/utils/types.h"

void print_matrix(char *file_name, const real_t *matrix, const int_t nrows,
    const int_t ncols);

void print_matrix_name(char *file_name, char *name, const real_t *matrix,
        const int_t nrows, const int_t ncols);

void print_array(char *file_name, real_t *array, int_t size);

#endif  // ACADOS_UTILS_PRINT_H_

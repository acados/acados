/*
 *    This file is part of ACADOS.
 *
 *    ACADOS is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADOS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADOS; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <stdio.h>
#include <string.h>
#include "acados/utils/print.h"

void print_matrix(char *file_name, const real_t *matrix, const int_t nrows,
    const int_t ncols) {

    FILE *output;
    if (strcmp(file_name, "stdout") == 0) {
        output = stdout;
    } else {
        output = fopen(file_name, "w");
    }
    if (output == NULL) {
        fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
    }
    // Assumes column major ordering
    for (int_t i = 0; i < nrows; i++) {
        for (int_t j = 0; j < ncols; j++) {
            fprintf(output, "%+.3e ", matrix[j*nrows+i]);
        }
        fprintf(output, "\n");
    }
    if (output != stdout) fclose(output);
}

void print_matrix_name(char *file_name, char *name, const real_t *matrix,
        const int_t nrows, const int_t ncols) {
    FILE *output;
    if (strcmp(file_name, "stdout") == 0) {
        output = stdout;
    } else {
        output = fopen(file_name, "w");
    }
    if (output == NULL) {
        fprintf(stderr, "Opening of file `%s' failed!\n", file_name);
    }
    fprintf(output, "%s:\n", name);
    // Assumes column major ordering
    for (int_t i = 0; i < nrows; i++) {
        for (int_t j = 0; j < ncols; j++) {
            fprintf(output, "%+.3e ", matrix[j*nrows+i]);
        }
        fprintf(output, "\n");
    }
    if (output != stdout) fclose(output);
}

void print_array(char *file_name, real_t *array, int_t size) {
    print_matrix(file_name, array, size, 1);
}

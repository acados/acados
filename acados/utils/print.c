/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acados/utils/print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
            fprintf(output, "%+.3e ", matrix[j * nrows + i]);
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
            fprintf(output, "%+.3e ", matrix[j * nrows + i]);
        }
        fprintf(output, "\n");
    }
    if (output != stdout) fclose(output);
}

void print_int_matrix(char *file_name, const int_t *matrix, const int_t nrows,
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
            fprintf(output, "%d ", matrix[j * nrows + i]);
        }
        fprintf(output, "\n");
    }
    if (output != stdout) fclose(output);
}

void print_array(char *file_name, real_t *array, int_t size) {
    print_matrix(file_name, array, size, 1);
}

void print_int_array(char *file_name, const int_t *array, int_t size) {
    print_int_matrix(file_name, array, size, 1);
}

// Read space delimited file into column-major matrix
void read_matrix(const char *file_name, real_t *array, const int_t nrows,
                 const int_t ncols) {
    FILE *file;
    file = fopen(file_name, "r");

    if (file == NULL) {
        printf("Error opening file %s ! ! ! ! ! ! ! ! !\n", file_name);
        exit(1);
    }

    // Read numbers from file into buffer.
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            if (!fscanf(file, "%lf", &array[nrows * j + i])) {
                break;
            }
        }
    }

    fclose(file);
}

void print_ocp_qp(ocp_qp_in *qp) {
    int_t N = qp->N;
    printf("ocp_qp structure with contents:\n");
    printf("N: %d\n", qp->N);
    printf("nx:\n");
    print_int_matrix("stdout", qp->nx, 1, N + 1);
    printf("nu:\n");
    print_int_matrix("stdout", qp->nu, 1, N + 1);
    printf("nb:\n");
    print_int_matrix("stdout", qp->nb, 1, N + 1);
    printf("nc:\n");
    print_int_matrix("stdout", qp->nc, 1, N + 1);
    for (int_t stage = 0; stage < N+1; stage++) {
        if (stage < N) {
            printf("A[%d]:\n", stage);
            print_matrix("stdout", qp->A[stage], qp->nx[stage], qp->nx[stage]);
            printf("B[%d]:\n", stage);
            print_matrix("stdout", qp->B[stage], qp->nx[stage], qp->nu[stage]);
            printf("b[%d]:\n", stage);
            print_matrix("stdout", qp->b[stage], qp->nx[stage], 1);
        }
        printf("Q[%d]:\n", stage);
        print_matrix("stdout", qp->Q[stage], qp->nx[stage], qp->nx[stage]);
        printf("R[%d]:\n", stage);
        print_matrix("stdout", qp->R[stage], qp->nu[stage], qp->nu[stage]);
        printf("S[%d]:\n", stage);
        print_matrix("stdout", qp->S[stage], qp->nu[stage], qp->nx[stage]);
        printf("q[%d]:\n", stage);
        print_matrix("stdout", qp->q[stage], qp->nx[stage], 1);
        printf("r[%d]:\n", stage);
        print_matrix("stdout", qp->r[stage], qp->nu[stage], 1);
        printf("lb[%d]:\n", stage);
        print_matrix("stdout", qp->lb[stage], qp->nb[stage], 1);
        printf("ub[%d]:\n", stage);
        print_matrix("stdout", qp->ub[stage], qp->nb[stage], 1);
    }
    printf("\n");
}

void print_ocp_qp_to_file(ocp_qp_in *qp) {
    char filename[MAX_STR_LEN];
    for (int_t i = 0; i <= qp->N; i++) {
        snprintf(filename, sizeof(filename), "Qm%d.txt", i);
        print_matrix(filename, qp->Q[i], qp->nx[i], qp->nx[i]);
        snprintf(filename, sizeof(filename), "Sm%d.txt", i);
        print_matrix(filename, qp->S[i], qp->nu[i], qp->nx[i]);
        snprintf(filename, sizeof(filename), "Rm%d.txt", i);
        print_matrix(filename, qp->R[i], qp->nu[i], qp->nu[i]);
        snprintf(filename, sizeof(filename), "qv%d.txt", i);
        print_matrix(filename, qp->q[i], qp->nx[i], 1);
        snprintf(filename, sizeof(filename), "rv%d.txt", i);
        print_matrix(filename, qp->r[i], qp->nu[i], 1);
        if (i < qp->N) {
            snprintf(filename, sizeof(filename), "Am%d.txt", i);
            print_matrix(filename, qp->A[i], qp->nx[i+1], qp->nx[i+1]);
            snprintf(filename, sizeof(filename), "Bm%d.txt", i);
            print_matrix(filename, qp->B[i], qp->nx[i+1], qp->nu[i]);
            snprintf(filename, sizeof(filename), "bv%d.txt", i);
            print_matrix(filename, qp->b[i], qp->nx[i+1], 1);
        }
        snprintf(filename, sizeof(filename), "idxb%d.txt", i);
        print_int_matrix(filename, qp->idxb[i], qp->nb[i], 1);
        snprintf(filename, sizeof(filename), "lb%d.txt", i);
        print_matrix(filename, qp->lb[i], qp->nb[i], 1);
        snprintf(filename, sizeof(filename), "ub%d.txt", i);
        print_matrix(filename, qp->ub[i], qp->nb[i], 1);
        snprintf(filename, sizeof(filename), "Cx%d.txt", i);
        print_matrix(filename, qp->Cx[i], qp->nc[i], qp->nx[i]);
        snprintf(filename, sizeof(filename), "Cu%d.txt", i);
        print_matrix(filename, qp->Cu[i], qp->nc[i], qp->nu[i]);
    }
}

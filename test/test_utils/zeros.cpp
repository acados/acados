#include "test/test_utils/zeros.h"
#include <cstdlib>

void i_zeros(int_t **pA, int_t row, int_t col) {
    void *temp = malloc((row*col)*sizeof(int_t));
    *pA = (int_t *) temp;
    int_t *A = *pA;
    int_t i;
    for (i = 0; i < row*col; i++) A[i] = 0;
}

void d_zeros(real_t **pA, int_t row, int_t col) {
    void *temp = malloc((row*col)*sizeof(real_t));
    *pA = (real_t *) temp;
    real_t *A = *pA;
    int_t i;
    for (i = 0; i < row*col; i++) A[i] = 0.0;
}

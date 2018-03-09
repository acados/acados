
#include <stdio.h>
#include <stdlib.h>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// make && gcc -g -o tst ../experimental/dimitris/main_ocp_qp.c -I../ -I../external -L ./acados -L ./external/blasfeo -lblasfeo -lacados && ./tst

#define N 10
#define NX 3
#define NU 2
#define NC 1

int main() {

    int_t nx[N+1];
    int_t nu[N+1];
    int_t nb[N+1];
    int_t nc[N+1];

    real_t A[NX*NX] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    for (int_t k = 0; k < N+1; k++) {
        nx[k] = NX;

        if (k < N) {
            nu[k] = NU;
        } else {
            nu[k] = 0;
        }

        nb[k] = nx[k] + nu[k];
        nc[k] = NC;
    }

    ocp_qp_in *qp_in = ocp_qp_in_create(N, nx, nu, nb, nc);

    printf("dimensions:\n");
    for (int_t k = 0; k < N+1; k++) {
        printf("nx[%d] = %d, nu[%d] = %d, nb[%d] = %d, nc[%d] = %d\n", k, qp_in->nx[k], k, qp_in->nu[k], k, qp_in->nb[k], k, qp_in->nc[k]);
    }

    real_t **hA = (real_t **) qp_in->A;

    for (int_t k = 0; k < N; k++)
        hA[k] = A;

    printf("\nA[%d] = \n", 3);
    d_print_mat(NX, NX, (real_t*) qp_in->A[3], NX);
    free(qp_in);
}
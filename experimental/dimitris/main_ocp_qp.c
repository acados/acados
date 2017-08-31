
#include <stdio.h>
#include <stdlib.h>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"

// make && gcc -g -o tst ../experimental/dimitris/main_ocp_qp.c -I../  -L ./acados -lacados && ./tst

int main() {

    int_t N = 10;
    int_t NX = 4;
    int_t NU = 2;
    int_t NC = 1;

    int_t nx[N+1];
    int_t nu[N+1];
    int_t nb[N+1];
    int_t nc[N+1];

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

    ocp_qp_in *qp_in = create_ocp_qp_in(N, nx, nu, nb, nc);

    printf("dimensions:\n");
    for (int_t k = 0; k < N+1; k++) {
        printf("nx[%d] = %d, nu[%d] = %d, nb[%d] = %d, nc[%d] = %d\n", k, qp_in->nx[k], k, qp_in->nu[k], k, qp_in->nb[k], k, qp_in->nc[k]);
    }

    free(qp_in);
}
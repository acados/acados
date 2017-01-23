#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_i_aux.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_d_blas.h"

#include "acados/sim/sim_collocation.h"
#include "acados/utils/print.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void get_Gauss_nodes(const int_t num_stages, real_t *nodes) {
//    if ( num_stages == 1 ) {         // GL2
//        nodes[0] = 1.0/2.0;
//    } else if ( num_stages == 2 ) {  // GL4
//        memcpy(nodes,
//                ((real_t[]) {1.0/2.0+sqrt(3.0)/6.0, 1.0/2.0-sqrt(3.0)/6.0}),
//                sizeof(*nodes) * (num_stages));
//    } else if ( num_stages == 3 ) {  // GL6
//        memcpy(nodes,
//                ((real_t[]) {1.0/2.0-sqrt(15.0)/10.0, 1.0/2.0, 1.0/2.0+sqrt(15.0)/10.0}),
//                sizeof(*nodes) * (num_stages));
//    } else {
//        // throw error somehow?
//    }
    uint N = num_stages-1;
    uint N1 = N+1;
    uint N2 = N+2;
    real_t *x_init;
    real_t *y;
    real_t *y_prev;
    real_t *lgvm;         // Legendre-Gauss Vandermonde Matrix
    real_t *der_lgvm;        // derivative of LGVM
    real_t err = 1;
    real_t eps = 2e-16;

    d_zeros(&x_init, N1, 1);
    d_zeros(&y, N1, 1);
    d_zeros(&y_prev, N1, 1);
    d_zeros(&lgvm, N1, N2);
    d_zeros(&der_lgvm, N1, 1);

    real_t a = 0.0;
    real_t b = 1.0;     // code for collocation interval [a,b]

    for (uint i = 0; i < N1; i++) {
        x_init[i] = -1+i*2.0/N;
        y[i] = cos((2*i+1)*M_PI/(2*N+2))+(0.27/N1)*sin(M_PI*x_init[i]*N/N2);
        y_prev[i] = 2.0;
    }

    while (err > eps) {     // iterate until step sufficiently small
        for (uint i = 0; i < N1; i++) lgvm[i] = 1.0;
        for (uint i = 0; i < N1; i++) lgvm[N1+i] = y[i];
        for (uint k = 2; k < N2; k++) {
            for (uint i = 0; i < N1; i++)
                lgvm[k*N1+i] = ((2*k-1)*y[i]*lgvm[(k-1)*N1+i]-(k-1)*lgvm[(k-2)*N1+i])/k;
        }
        for (uint i = 0; i < N1; i++)
            der_lgvm[i] = N2*(lgvm[N*N1+i]-y[i]*lgvm[N1*N1+i])/(1-pow(y[i], 2));
        for (uint i = 0; i < N1; i++) y_prev[i] = y[i];

        for (uint i = 0; i < N1; i++)
            y[i] = y_prev[i] - lgvm[N1*N1+i]/der_lgvm[i];     // Newton step
//        print_matrix("stdout", y, 1, num_stages);

        err = 0;
        for (uint i = 0; i < N1; i++) {
            if (err < fabs(y[i]-y_prev[i])) err = fabs(y[i]-y_prev[i]);
        }
    }
    for (uint i = 0; i < N1; i++) nodes[i] = (a*(1-y[i]) + b*(1+y[i]))/2;

//    free(x_init);
//    free(y);
//    free(y_prev);
//    free(lgvm);
//    free(der_lgvm);
}


void create_Butcher_table(const int_t num_stages, const real_t *nodes,
        real_t *b, real_t *A) {
        b[0] = nodes[0];
//    if ( strcmp(name, "Gauss") == 0 ) {  // GAUSS METHODS
        if ( num_stages == 1 ) {         // GL2
            A[0] = 1.0/2.0;
            b[0] = 1.0;
        } else if ( num_stages == 2 ) {  // GL4
            memcpy(A,
                    ((real_t[]) {1.0/4.0, (1.0/4.0-sqrt(3.0)/6.0),
                    (1.0/4.0+sqrt(3.0)/6.0), 1.0/4.0}),
                    sizeof(*A) * (num_stages*num_stages));
            memcpy(b,
                    ((real_t[]) {1.0/2.0, 1.0/2.0}), sizeof(*b) * (num_stages));
        } else if ( num_stages == 3 ) {  // GL6
            memcpy(A,
                    ((real_t[]) {5.0/36.0, 5.0/36.0+1.0/24.0*sqrt(15.0),
                    5.0/36.0+1.0/30.0*sqrt(15.0), 2.0/9.0-1.0/15.0*sqrt(15.0),
                    2.0/9.0, 2.0/9.0+1.0/15.0*sqrt(15.0), 5.0/36.0-1.0/30.0*sqrt(15.0),
                    5.0/36.0-1.0/24.0*sqrt(15.0), 5.0/36.0}),
                    sizeof(*A) * (num_stages*num_stages));
            memcpy(b,
                    ((real_t[]) {5.0/18.0, 4.0/9.0, 5.0/18.0}),
                    sizeof(*b) * (num_stages));
        } else {
            // throw error somehow?
        }
//    } else if ( strcmp(name, "Radau") == 0 ) {
//        // TODO(rien): add Radau IIA collocation schemes
//    } else {
//        // throw error somehow?
//    }
}

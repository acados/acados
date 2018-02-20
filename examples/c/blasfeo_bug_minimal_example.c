#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// blasfeo
#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>
int main() {
    int nx = 8;
    int nu = 2;
    int nZ = 4;
    int nK1 = nx * 4;
    struct blasfeo_dmat A;
    struct blasfeo_dmat B;
    struct blasfeo_dmat C;
    struct blasfeo_dmat result;

    double *some_doubles;
    some_doubles = (double*) calloc(10, sizeof(double));
    for (int ii = 0; ii < 10; ii++) {
        some_doubles[ii] = (double) ii;
    }
    blasfeo_allocate_dmat(nx, nx+nu, &B); //    
    blasfeo_allocate_dmat(nK1, nx, &A); //    
    blasfeo_allocate_dmat(nK1, nu, &result); //
    blasfeo_allocate_dmat(nK1, nu, &C); //   

    for (int ii = 0; ii < nx; ii++) {
        blasfeo_pack_dmat(1,10, &some_doubles[0], 1, &B, ii,0);
    }
    blasfeo_dgese(nK1, nx, 1.0, &A, 0, 0);

    printf("A = \n");
    blasfeo_print_dmat(nK1,  nx, &A, 0, 0);
    printf("B_multiplication = \n");
    blasfeo_print_dmat(nx, nu, &B, 0, nx);    
    blasfeo_dgemm_nn(nK1, nu,  nx, -1.0, &A, 0, 0, &B, 0, nx, 1.0, &C, 0, 0, &result , 0, 0); // Blasfeo HP & Reference differ here

    printf("A * B_multiplication: result = \n");
    blasfeo_print_exp_dmat(nK1,  nu, &result,0,0);

    return 0;
}
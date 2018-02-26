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
    int nx = 9;
    int nu = 2;

    struct blasfeo_dmat A;
    struct blasfeo_dvec lambda;
    struct blasfeo_dvec lambda0;
    struct blasfeo_dvec lambda_new;

    blasfeo_allocate_dmat(nx, nx, &A); //    
    blasfeo_allocate_dvec(nx+nu, &lambda);
    blasfeo_allocate_dvec(nx+nu, &lambda0);
    blasfeo_allocate_dvec(nx+nu, &lambda_new);
    for (int ii = 0; ii < nx+nu; ii++) {
        blasfeo_dvecin1((double) ii, &lambda, ii);
    }
    blasfeo_dgese(nx, nx, 1.0, &A, 0, 0);

    blasfeo_print_dmat(nx, nx, &A, 0, 0);
    printf("lambda = \n");
    blasfeo_print_dvec(nx+nu, &lambda, 0);    

    blasfeo_dgemv_t(nx, nx, 1.0, &A, 0, 0, &lambda, 0, 0.0, &lambda0, 0, &lambda_new, 0); // recheck!
    printf("lambda: result = \n");
    blasfeo_print_exp_dvec(nx+  nu, &lambda_new,0);

    return 0;
}
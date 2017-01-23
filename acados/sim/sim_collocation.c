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

void get_Gauss_nodes(const int_t num_stages, real_t *nodes) {
    if ( num_stages == 1 ) {         // GL2
        nodes[0] = 1.0/2.0;
    } else if ( num_stages == 2 ) {  // GL4
        memcpy(nodes,
                ((real_t[]) {1.0/2.0+sqrt(3.0)/6.0, 1.0/2.0-sqrt(3.0)/6.0}),
                sizeof(*nodes) * (num_stages));
    } else if ( num_stages == 3 ) {  // GL6
        memcpy(nodes,
                ((real_t[]) {1.0/2.0-sqrt(15.0)/10.0, 1.0/2.0, 1.0/2.0+sqrt(15.0)/10.0}),
                sizeof(*nodes) * (num_stages));
    } else {
        // throw error somehow?
    }
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

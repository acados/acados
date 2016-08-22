#ifndef ACADOS_ACADOS_TYPES_H_
#define ACADOS_ACADOS_TYPES_H_

typedef double real_t;
typedef unsigned int uint;
typedef int int_t;

#define NNN 20
#define NX 8
#define NU 3

// enum of return values
enum return_values{
    ACADOS_SUCCESS,
    ACADOS_MAXITER,
    ACADOS_MINSTEP
};

#endif  // ACADOS_ACADOS_TYPES_H_

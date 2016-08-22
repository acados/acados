#ifndef ACADOS_ACADOS_TYPES_H_
#define ACADOS_ACADOS_TYPES_H_

typedef double real_t;
typedef unsigned int uint;
typedef int int_t;

#define NNN 13
#define NX 2
#define NU 1

// enum of return values
enum return_values{
    ACADOS_SUCCESS,
    ACADOS_MAXITER,
    ACADOS_MINSTEP
};

#endif  // ACADOS_ACADOS_TYPES_H_

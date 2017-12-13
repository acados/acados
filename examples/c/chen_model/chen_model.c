/* The model comes from \cite{Chen1998} */
#include "examples/c/chen_model/chen_model.h"

/* Ignore vde */
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
#if __GNUC__ >= 6
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#else
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#endif

#define LINEAR_MODEL 0
#define NX 2
#define NU 1

void VDE_fun(const int_t nx, const int_t nu, const real_t* in, real_t* out,
             int (*vde)(const real_t**, real_t**, int*, real_t*, int)) {
    const real_t* x = in;
    const real_t* u = in + NX + NX * (NX + NU);

#if LINEAR_MODEL == 0
    /* COMPUTE AUXILIARY VARIABLES: */
    /* ---------------------------- */
    real_t aux[12];
    aux[0] = ((u[0] * (real_t)(0.5)) * x[2]);
    aux[1] = (aux[0] + x[4]);
    aux[2] = ((u[0] * (real_t)(0.5)) * x[3]);
    aux[3] = (aux[2] + x[5]);
    aux[4] = x[2];
    aux[5] = (aux[4] + ((u[0] * ((real_t)(0.) - (real_t)(2.))) * x[4]));
    aux[6] = x[3];
    aux[7] = (aux[6] + ((u[0] * ((real_t)(0.) - (real_t)(2.))) * x[5]));
    aux[8] = ((u[0] * (real_t)(0.5)) * x[6]);
    aux[9] = (aux[8] + x[7]);
    aux[10] = x[6];
    aux[11] = (aux[10] + ((u[0] * ((real_t)(0.) - (real_t)(2.))) * x[7]));

    /* COMPUTE OUTPUTS: */
    /* ---------------- */
    out[0] = (x[1] + (u[0] * ((real_t)(0.5) + ((real_t)(0.5) * x[0]))));
    out[1] = (x[0] + (u[0] * ((real_t)(0.5) - ((real_t)(2.) * x[1]))));
    out[2] = aux[1];
    out[3] = aux[3];
    out[4] = aux[5];
    out[5] = aux[7];
    out[6] = (aux[9] + ((real_t)(0.5) + ((real_t)(0.5) * x[0])));
    out[7] = (aux[11] + ((real_t)(0.5) - ((real_t)(2.) * x[1])));
#else
    /* COMPUTE OUTPUTS: */
    /* ---------------- */
    out[0] = x[1];
    out[1] = u[0];
    out[2] = x[4];
    out[3] = x[5];
    out[4] = 0;
    out[5] = 0;
    out[6] = x[7];
    out[7] = 1;
#endif
}

void jac_fun(const real_t* in, real_t* out) {
    const real_t* x = in;
    const real_t* u = in + NX;

    /* Compute outputs: */
    out[0] = (u[0] * (real_t)(0.5));
    out[1] = (real_t)(1.0);
    out[2] = ((real_t)(0.5) + ((real_t)(0.5) * x[0]));
    out[3] = (real_t)(1.0);
    out[4] = (u[0] * (-(real_t)(2.0)));
    out[5] = ((real_t)(0.5) - ((real_t)(2.0) * x[1]));
}

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#if __GNUC__ >= 6
#pragma GCC diagnostic pop
#endif
#endif

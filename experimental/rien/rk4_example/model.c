/* The model comes from \cite{Chen1998} */

#include "common_header.h"

real_t aux[12];

void VDE_fun( const real_t* in, real_t* out ){
const real_t* x = in;
const real_t* u  = in + NX + NX*(NX+NU);

/* COMPUTE AUXILIARY VARIABLES: */
/* ---------------------------- */
aux[0] = ((u[0]*(real_t)(0.5))*x[2]);
aux[1] = (aux[0]+x[4]);
aux[2] = ((u[0]*(real_t)(0.5))*x[3]);
aux[3] = (aux[2]+x[5]);
aux[4] = x[2];
aux[5] = (aux[4]+((u[0]*((real_t)(0.)-(real_t)(2.)))*x[4]));
aux[6] = x[3];
aux[7] = (aux[6]+((u[0]*((real_t)(0.)-(real_t)(2.)))*x[5]));
aux[8] = ((u[0]*(real_t)(0.5))*x[6]);
aux[9] = (aux[8]+x[7]);
aux[10] = x[6];
aux[11] = (aux[10]+((u[0]*((real_t)(0.)-(real_t)(2.)))*x[7]));

/* COMPUTE OUTPUTS: */
/* ---------------- */
out[0] = (x[1]+(u[0]*((real_t)(0.5)+((real_t)(0.5)*x[0]))));
out[1] = (x[0]+(u[0]*((real_t)(0.5)-((real_t)(2.)*x[1]))));
out[2] = aux[1];
out[3] = aux[3];
out[4] = aux[5];
out[5] = aux[7];
out[6] = (aux[9]+((real_t)(0.5)+((real_t)(0.5)*x[0])));
out[7] = (aux[11]+((real_t)(0.5)-((real_t)(2.)*x[1])));
}

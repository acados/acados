#ifndef AIRCRAFT_INTEGRATOR_H_
#define AIRCRAFT_INTEGRATOR_H_

#include "acados/utils/types.h"

typedef struct _rk4_int{

  real_t* x_in;
  real_t* u_in;
  real_t* p_in;
  real_t t0_in;
  real_t h_in;
  real_t N_in;

  real_t nx;
  real_t nu;
  real_t np;
  real_t n_steps;

  real_t *x_out;
  real_t *Sx_out;
  real_t *Su_out;

  void (*eval_dynamics)(const real_t**, real_t**, int*, real_t*, int) ;
  void (*eval_dynamics_work)(int *, int* , int *, int *);

  // additional internal buffers
  real_t *intermediate_Sx_in;
  real_t *intermediate_Sx_out;
  real_t *intermediate_x_in;
  real_t *intermediate_x_out;
  real_t *intermediate_Su;
  real_t *temp_Sx_out;
  real_t *temp_Su_out;

  real_t *internal_mem;

} rk4_int;

void integrate_aircraft_ode(rk4_int* rk4_int);

#endif  // AIRCRAFT_INTEGRATOR_H_

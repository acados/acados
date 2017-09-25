/* The model comes from \cite{Wirsching2006} */
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"

#include "../models/aircraft_integrator.h"
#include "../acados/utils/tools.h"

// The auto-generated VDE functions from CasADi:
void integrate_aircraft_ode(rk4_int* rk4_int) {

    const int_t NX = rk4_int->nx;
    const int_t NU = rk4_int->nu;
    const int_t NP = rk4_int->np;

    real_t* x_out = rk4_int->x_out;
    real_t* Sx_out = rk4_int->Sx_out;
    real_t* Su_out = rk4_int->Su_out;

    int casadi_mem = 0;
    int *casadi_iw = 0;
    real_t *casadi_w = rk4_int->internal_mem;

    const real_t *casadi_arg[5];
    real_t *casadi_res[3];

    real_t *intermediate_x_in = rk4_int->intermediate_x_in;
    real_t *intermediate_x_out = rk4_int->intermediate_x_out;

    for (int_t i = 0; i < NX; i++) {
      intermediate_x_in[i] = rk4_int->x_in[i];
      intermediate_x_out[i] = rk4_int->x_in[i];
    }

    casadi_arg[0] = intermediate_x_in;
    casadi_arg[1] = rk4_int->u_in;
    casadi_arg[2] = rk4_int->p_in;
    casadi_arg[3] = &rk4_int->t0_in;
    real_t h = rk4_int->h_in/rk4_int->n_steps;
    casadi_arg[4] = &h;

    // stuff should be initialized already, but never hurts...
    real_t *intermediate_Sx_in = rk4_int->intermediate_Sx_in;
    for (int_t i = 0; i < NX*NX; i++) intermediate_Sx_in[i] = 0.0;

    real_t *intermediate_Sx_out = rk4_int->intermediate_Sx_out;
    for (int_t i = 0; i < NX*NX; i++) intermediate_Sx_out[i] = 0.0;

    real_t *intermediate_Su = rk4_int->intermediate_Su;
    for (int_t i = 0; i < NX*NU; i++) intermediate_Su[i] = 0.0;

    real_t *temp_Sx_out = rk4_int->temp_Sx_out;
    for (int_t i = 0; i < NX*NX; i++) temp_Sx_out[i] = 0.0;

    real_t *temp_Su_out = rk4_int->temp_Su_out;
    for (int_t i = 0; i < NX*NU; i++) temp_Su_out[i] = 0.0;

    for (int_t i = 0; i < NX; i++) intermediate_Sx_in[NX*i + i] = 1.0;

    casadi_res[2] = intermediate_x_out;  // Andrea: sensitivies first and then states
    casadi_res[0] = temp_Sx_out;
    casadi_res[1] = temp_Su_out;

    // copy intermediate_x_out
    for (int_t i = 0; i < NX; i++) x_out[i] = intermediate_x_out[i];

    for (int_t i = 0; i < rk4_int->n_steps; i++) {

      rk4_int->eval_dynamics(casadi_arg, casadi_res, casadi_iw, casadi_w, casadi_mem);

      // propagate sensitivities
      dgemm_nn_3l(NX, NX, NX, temp_Sx_out, NX, intermediate_Sx_in, NX, intermediate_Sx_out, NX);
      dgemm_nn_3l(NX, NU, NX, temp_Sx_out, NX, intermediate_Su, NX, temp_Su_out, NX);

      // copy result to intermediate_Sx_in and intermediate_Sx_out
      for (int_t j = 0; j < NX*NX; j++) {
        intermediate_Sx_in[j] = intermediate_Sx_out[j];
        intermediate_Sx_out[j] = 0.0;
      }

      for (int_t j = 0; j < NX*NU; j++) intermediate_Su[j] = temp_Su_out[j];

      // update state
      for (int_t i = 0; i < NX; i++) intermediate_x_in[i] = intermediate_x_out[i];
    }

    // copy results to output buffers
    // copy intermediate_Sx_in and intermediate_Sx_out
    for (int_t i = 0; i < NX*NX; i++) Sx_out[i] = intermediate_Sx_in[i];
    for (int_t i = 0; i < NX*NU; i++) Su_out[i] = intermediate_Su[i];

    // copy intermediate_x_out
    for (int_t i = 0; i < NX; i++) x_out[i] = intermediate_x_out[i];

}

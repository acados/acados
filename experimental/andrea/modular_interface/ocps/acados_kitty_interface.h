#ifndef KITTY_ACADOS_INTERFACE_H_
#define KITTY_ACADOS_INTERFACE_H_

// ACADOS headers
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/utils/tools.h"

#include "acados/sim/sim_erk_integrator.h"
#include "../models/aircraft_integrator.h"

#include <stdio.h>

typedef struct _init_data {  // data needed to allocate memory (cannot be changed online)

  int_t NN;  // number  of stages
  int_t MM;  // number  of stages
  int_t NX;  // number  of states
  int_t NU;  // number  of inputs
  int_t NB0; // number  of b  ounds on stage 0
  int_t NB;  // number  of bounds on stage 1 to N-1
  int_t NBN; // number  of bounds on stage N
  int_t NP;  // number  of parameters
  int_t NR;  // number  of residuals (least-square)

  int_t *idxb0; //bounds index on stage 0
  int_t *idxb;  //bounds index on stage 1 to N-1
  int_t *idxbN; //bounds index on stage N

  double sigma_mu;

} init_nmpc_data;

typedef struct _nmpc_data {

  // ocp nlp data
  real_t *w;
  real_t *lam_n;
  real_t *t_n;

  int_t max_sqp_iters;
  int_t max_iters;

  int_t NN;  // number  of stages
  int_t MM;  // number  of stages (dual MM)
  int_t NX;  // number  of states
  int_t NU;  // number  of inputs
  int_t NB0; // number  of bounds on stage 0
  int_t NB;  // number  of bounds on stage 1 to N-1
  int_t NBN; // number  of bounds on stage N
  int_t NP;  // number  of parameters
  int_t NR;  // number  of residuals (least-square)

  int_t *idxb0; //bounds index on stage 0
  int_t *idxb;  //bounds index on stage 1 to N-1
  int_t *idxbN; //bounds index on stage N


  real_t *lb0; //bounds on stage 0
  real_t *ub0;

  real_t *lb; //bounds on stage N-1
  real_t *ub;

  real_t *lbN; //bounds on stage N
  real_t *ubN;

  real_t regQ;
  real_t regR;

  void (*res_x)(const real_t**, real_t**, int*, real_t*, int);  // ocp_xtracing_something
  void (*res_x_work)(int *, int* , int *, int *);               // ocp_xtracing_something
  void (*res_u)(const real_t**, real_t**, int*, real_t*, int);  // ocp_utracing_something
  void (*res_u_work)(int *, int* , int *, int *);               // ocp_utracing_something

  // stuff you don't want to care about
  void *residual_x_eval_mem;
  real_t *residual_x_out;
  real_t *residual_x_in;
  real_t *drdx_tran;

  void *residual_u_eval_mem;
  real_t *residual_u_out;
  real_t *residual_u_in;
  real_t *drdu_tran;

  // qp data
  ocp_qp_in *qp_in;
  ocp_qp_out *qp_out;

  // solver data
  ocp_qp_hpmpc_args *hpmpc_args;
  void *workspace;

  // gnuplot stuff
  FILE *gnuplotPipe;

} nmpc_data;


typedef struct _acados_options {

  int print_level;
  real_t non_symmetry_tol;
  real_t newton_step_size;
  int_t sqp_steps;
  int_t plot_open_loop;
  int_t shifting;
  real_t max_qp_step;
  real_t terminal_cost_scaling;

  // linear vs nonlinear residuals
  int nls;

  int_t use_gnuplot;

} acados_options;

#if defined (__cplusplus)
extern "C" {
#endif

int_t run_acados(nmpc_data *nmpc_data, rk4_int* rk4_int, acados_options* acados_options);
void init_acados(nmpc_data *nmpc_data, rk4_int* rk4_int, init_nmpc_data* init_data, acados_options* acados_options);
int_t de_init_acados(nmpc_data *nmpc_data, rk4_int* rk4_int);

#if defined (__cplusplus)
}
#endif

#endif   // ACADOS_INTERFACE_H_

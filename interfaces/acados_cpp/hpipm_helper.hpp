
#ifndef ACADOS_INTERFACES_ACADOS_CPP_HPIPM_HELPER_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_HPIPM_HELPER_HPP_

#define num_rows_Q(stage, dim) (dim->nx[stage])

#define num_cols_Q(stage, dim) (dim->nx[stage])

#define num_rows_S(stage, dim) (dim->nu[stage])

#define num_cols_S(stage, dim) (dim->nx[stage])

#define num_rows_R(stage, dim) (dim->nu[stage])

#define num_cols_R(stage, dim) (dim->nu[stage])

#define num_elems_q(stage, dim) (dim->nx[stage])

#define num_elems_r(stage, dim) (dim->nu[stage])

#define num_rows_A(stage, dim) (dim->nx[stage+1])

#define num_cols_A(stage, dim) (dim->nx[stage])

#define num_rows_B(stage, dim) (dim->nx[stage+1])

#define num_cols_B(stage, dim) (dim->nu[stage])

#define num_elems_b(stage, dim) (dim->nx[stage+1])

#define num_elems_lbx(stage, dim) (dim->nbx[stage])

#define num_elems_lbu(stage, dim) (dim->nbu[stage])

#define num_elems_ubx(stage, dim) (dim->nbx[stage])

#define num_elems_ubu(stage, dim) (dim->nbu[stage])

#define num_rows_C(stage, dim) (dim->ng[stage])

#define num_cols_C(stage, dim) (dim->nx[stage])

#define num_rows_D(stage, dim) (dim->ng[stage])

#define num_cols_D(stage, dim) (dim->nu[stage])

#define num_elems_lg(stage, dim) (dim->ng[stage])

#define num_elems_ug(stage, dim) (dim->ng[stage])

#endif  // ACADOS_INTERFACES_ACADOS_CPP_HPIPM_HELPER_HPP_

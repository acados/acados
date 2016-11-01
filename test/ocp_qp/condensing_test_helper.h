#ifndef TEST_OCP_QP_CONDENSING_TEST_HELPER_H_
#define TEST_OCP_QP_CONDENSING_TEST_HELPER_H_

#include "acados/ocp_qp/condensing.h"

void calculate_num_state_bounds(const ocp_qp_in *in, condensing_workspace *work);

int_t get_num_condensed_vars(const ocp_qp_in *in);

int_t get_num_constraints(const ocp_qp_in *in, condensing_workspace *work);

void allocateForUnconstrainedQPData(int_t N, int_t nx, int_t nu, ocp_qp_in * const qp);

void allocateForCondensingData(const ocp_qp_in * const qp_in, condensing_in *in,
    condensing_out *out, condensing_workspace *work);

#endif  // TEST_OCP_QP_CONDENSING_TEST_HELPER_H_

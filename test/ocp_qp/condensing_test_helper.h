#ifndef TEST_OCP_QP_CONDENSING_TEST_HELPER_H_
#define TEST_OCP_QP_CONDENSING_TEST_HELPER_H_

#include <string>
#include "test/test_utils/eigen.h"
#include "acados/ocp_qp/condensing.h"

void calculate_num_state_bounds(const ocp_qp_in *in, condensing_workspace *work);

int_t get_num_condensed_vars(const ocp_qp_in *in);

int_t get_num_constraints(const ocp_qp_in *in, condensing_workspace *work);

void allocateUnconstrainedQPData(int_t N, int_t nx, int_t nu, ocp_qp_in * const qp);

void allocateCondensingData(const ocp_qp_in * const qp_in, condensing_in *in,
    condensing_out *out, condensing_workspace *work);

void fillWithUnconstrainedData(ocp_qp_in *qp, Eigen::VectorXd *x0, std::string scenario);

void fillWithBoundsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu, std::string scenario);

void fillWithGeneralConstraintsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu,
    std::string scenario);

#endif  // TEST_OCP_QP_CONDENSING_TEST_HELPER_H_

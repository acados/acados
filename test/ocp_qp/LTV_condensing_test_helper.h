#ifndef TEST_OCP_QP_LTV_CONDENSING_TEST_HELPER_H_
#define TEST_OCP_QP_LTV_CONDENSING_TEST_HELPER_H_

#include "test/test_utils/eigen.h"
#include "acados/ocp_qp/condensing.h"

void readTVUnconstrainedInputDataFromFile(int_t nx, int_t nu,
    Eigen::MatrixXd *A, Eigen::MatrixXd *B, Eigen::VectorXd *b, Eigen::VectorXd *x0,
    Eigen::MatrixXd *Q, Eigen::MatrixXd *S, Eigen::MatrixXd *R,
    Eigen::VectorXd *q, Eigen::VectorXd *r);

void fillWithTVUnconstrainedData(ocp_qp_in *qp, Eigen::VectorXd *x0);

void readTVBoundsDataFromFile(int_t nx, int_t nu, Eigen::VectorXd *lbwx,
    Eigen::VectorXd *ubwx, Eigen::VectorXd *lbwu, Eigen::VectorXd *ubwu);

void fillWithTVBoundsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu);

void readTVGeneralConstraintsDataFromFile(int_t nx, int_t nu, int_t nc,
    Eigen::MatrixXd *Cx, Eigen::MatrixXd *Cu, Eigen::VectorXd *lbc, Eigen::VectorXd *ubc);

void fillWithTVGeneralConstraintsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu);

#endif  // TEST_OCP_QP_LTV_CONDENSING_TEST_HELPER_H_

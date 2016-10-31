#include "Eigen/Dense"
#include "acados/ocp_qp/condensing.h"

void readUnconstrainedInputDataFromFile(int_t nx, int_t nu,
    Eigen::MatrixXd *A, Eigen::MatrixXd *B, Eigen::VectorXd *b, Eigen::VectorXd *x0,
    Eigen::MatrixXd *Q, Eigen::MatrixXd *S, Eigen::MatrixXd *R,
    Eigen::VectorXd *q, Eigen::VectorXd *r);

void fillWithUnconstrainedData(ocp_qp_in *qp, Eigen::VectorXd *x0);

void readBoundsDataFromFile(int_t nx, int_t nu, Eigen::VectorXd *lbwx,
    Eigen::VectorXd *ubwx, Eigen::VectorXd *lbwu, Eigen::VectorXd *ubwu);

void fillWithBoundsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu);

void readGeneralConstraintsDataFromFile(int_t nx, int_t nu, int_t nc,
    Eigen::MatrixXd *Cx, Eigen::MatrixXd *Cu, Eigen::VectorXd *lbc, Eigen::VectorXd *ubc);

void fillWithGeneralConstraintsData(ocp_qp_in *qp, int_t N, int_t nx, int_t nu);

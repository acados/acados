#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>

#include "acados_cpp/ocp_qp.h"

#include "acados_c/ocp_qp.h"
#include "acados/utils/print.h"
#include "blasfeo/include/blasfeo_d_aux.h"


OcpQp::OcpQp(int N, std::vector<int> nx, std::vector<int> nu, std::vector<int> nb,
             std::vector<int> nc) : N(N), dimensions(nullptr), qp(nullptr) {

    if (N <= 0) throw std::invalid_argument("Number of stages must be positive");

    dimensions = create_ocp_qp_dims(N);
    dimensions->N = N;

    copyDimensions(nx, dimensions->nx);
    copyDimensions(nu, dimensions->nu);
    copyDimensions(nb, dimensions->nb);
    copyDimensions(nc, dimensions->ng);

    qp = create_ocp_qp_in(dimensions);
}

std::ostream& operator<<(std::ostream& oss, const OcpQp& qp) {
    print_ocp_qp_in(qp.qp);
    return oss;
}

void OcpQp::setQ(int i, double *Q) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_rows = dimensions->nx[i], num_cols = dimensions->nx[i];
    int row_offset = dimensions->nu[i], col_offset = dimensions->nu[i];
    blasfeo_pack_dmat(num_rows, num_cols, Q, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
}

void OcpQp::setS(int i, double *S) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    // In HPIPM, S is an nu times nx matrix
    int num_rows = dimensions->nu[i], num_cols = dimensions->nx[i];
    int row_offset = dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, S, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
}

void OcpQp::setR(int i, double *R) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_rows = dimensions->nu[i], num_cols = dimensions->nu[i];
    int row_offset = 0, col_offset = 0;
    blasfeo_pack_dmat(num_rows, num_cols, R, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
}

void OcpQp::setr(int i, double *r) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_rows = dimensions->nu[i], num_cols = 1;
    int row_offset = dimensions->nu[i] + dimensions->nx[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, r, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
    blasfeo_pack_dvec(num_rows, r, &(qp->rq[i]), 0);
}

void OcpQp::setq(int i, double *q) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_rows = dimensions->nx[i], num_cols = 1;
    int row_offset = dimensions->nu[i] + dimensions->nx[i], col_offset = dimensions->nu[i];
    blasfeo_pack_tran_dmat(num_rows, num_cols, q, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
    blasfeo_pack_dvec(num_rows, q, &(qp->rq[i]), col_offset);
}

void OcpQp::setA(int i, double *A) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + ").");
    int num_rows = dimensions->nx[i+1], num_cols = dimensions->nx[i];
    int row_offset = dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, A, num_rows, &(qp->BAbt[i]), row_offset, col_offset);
}

void OcpQp::setB(int i, double *B) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + ").");
    int num_rows = dimensions->nx[i+1], num_cols = dimensions->nu[i];
    int row_offset = 0, col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, B, num_rows, &(qp->BAbt[i]), row_offset, col_offset);
}

void OcpQp::setb(int i, double *b) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + ").");
    int num_rows = dimensions->nx[i+1], num_cols = 1;
    int row_offset = dimensions->nx[i] + dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, b, num_rows, &(qp->BAbt[i]), row_offset, col_offset);
    blasfeo_pack_dvec(num_rows, b, &(qp->b[i]), 0);
}

void OcpQp::setlb(int i, double *lb) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_elems = dimensions->nb[i];
    int offset = 0;
    blasfeo_pack_dvec(num_elems, lb, &(qp->d[i]), offset);
}

void OcpQp::setub(int i, double *ub) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_elems = dimensions->nb[i];
    int offset = dimensions->nb[i] + dimensions->ng[i];
    blasfeo_pack_dvec(num_elems, ub, &(qp->d[i]), offset);
    blasfeo_dvecsc(num_elems, -1.0, &(qp->d[i]), offset);
}

void OcpQp::setC(int i, double *C) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_rows = dimensions->nx[i], num_cols = dimensions->ng[i];
    int row_offset = dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, C, num_rows, &(qp->DCt[i]), row_offset, col_offset);
}

void OcpQp::setD(int i, double *D) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_rows = dimensions->nu[i], num_cols = dimensions->ng[i];
    int row_offset = 0, col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, D, num_rows, &(qp->DCt[i]), row_offset, col_offset);
}

void OcpQp::setlg(int i, double *lg) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_elems = dimensions->ng[i];
    int offset = dimensions->nb[i];
    blasfeo_pack_dvec(num_elems, lg, &(qp->d[i]), offset);
}

void OcpQp::setug(int i, double *ug) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" + std::to_string(N) + "].");
    int num_elems = dimensions->ng[i];
    int offset = 2*dimensions->nb[i] + dimensions->ng[i];
    blasfeo_pack_dvec(num_elems, ug, &(qp->d[i]), offset);
    blasfeo_dvecsc(num_elems, -1.0, &(qp->d[i]), offset);
}

void OcpQp::copyDimensions(std::vector<int> dimensions, int *dimension_ptr) {
    if (dimensions.size() != N+1)
        throw std::invalid_argument("Dimensions must be defined for all N+1 nodes.");
    if (std::any_of(dimensions.begin(), dimensions.end(), [](int a) {return a < 0;}))
        throw std::invalid_argument("Dimension must be non-negative for all N+1 nodes.");
    std::copy(dimensions.begin(), dimensions.end(), dimension_ptr);
}


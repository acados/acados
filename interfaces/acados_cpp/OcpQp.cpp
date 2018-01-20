#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <stdexcept>

#include "acados_cpp/OcpQp.h"

#include "acados_c/ocp_qp.h"
#include "acados/utils/print.h"
#include "blasfeo/include/blasfeo_d_aux.h"


OcpQp::OcpQp(int N, std::vector<int> nx, std::vector<int> nu, std::vector<int> nbx, 
             std::vector<int> nbu, std::vector<int> nc) : N(N), dimensions(nullptr), qp(nullptr) {

    if (N <= 0) throw std::invalid_argument("Number of stages must be positive");

    dimensions = create_ocp_qp_dims(N);
    dimensions->N = N;

    copyDimensions(nx, dimensions->nx);
    copyDimensions(nu, dimensions->nu);
    copyDimensions(nbx, dimensions->nbx);
    copyDimensions(nbu, dimensions->nbu);
    std::vector<int> nb(N+1);
    for(int i = 0; i <= N; i++)
        nb.at(i) = nbx.at(i) + nbu.at(i);
    copyDimensions(nb, dimensions->nb);
    copyDimensions(nc, dimensions->ng);
    std::vector<int> ns(N+1);
    copyDimensions(ns, dimensions->ns);

    qp = create_ocp_qp_in(dimensions);
}

OcpQp::OcpQp(int N, int nx, int nu, int nbx, int nbu, int nc)
    : OcpQp(N, std::vector<int>(N+1, nx), std::vector<int>(N+1, nu), std::vector<int>(N+1, nbx),
      std::vector<int>(N+1, nbu), std::vector<int>(N+1, nc)) {}

std::ostream& operator<<(std::ostream& oss, const OcpQp& qp) {
    print_ocp_qp_in(qp.qp);
    return oss;
}

std::vector<double> OcpQp::getA(int i) {
    int num_rows = dimensions->nx[i+1], num_cols = dimensions->nx[i];
    double A_col_major[num_rows * num_cols];
    blasfeo_unpack_tran_dmat(num_rows, num_cols, &(qp->BAbt[i]), dimensions->nu[i], 0,
                             A_col_major, num_rows);
    std::vector<double> A_raveled(num_rows * num_cols);
    std::copy_n(A_col_major, num_rows * num_cols, A_raveled.begin());
    return A_raveled;
}

int OcpQp::numRowsA(int i) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    return dimensions->nx[i+1];
}

int OcpQp::numColsA(int i) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    return dimensions->nx[i];
}

int OcpQp::numRowsB(int i) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    return dimensions->nx[i+1];
}

int OcpQp::numColsB(int i) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    return dimensions->nu[i];
}

int OcpQp::numRowsb(int i) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    return dimensions->nx[i+1];
}

int OcpQp::numColsb(int i) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    return 1;
}

int OcpQp::numRowsQ(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nx[i];
}

int OcpQp::numColsQ(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nx[i];
}

int OcpQp::numRowsS(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nu[i];
}

int OcpQp::numColsS(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nx[i];
}

int OcpQp::numRowsR(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nu[i];
}

int OcpQp::numColsR(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nu[i];
}

int OcpQp::numRowsq(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nx[i];
}

int OcpQp::numColsq(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return 1;
}

int OcpQp::numRowsr(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nu[i];
}

int OcpQp::numColsr(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return 1;
}

int OcpQp::numRowslb(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nb[i];
}

int OcpQp::numColslb(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return 1;
}

int OcpQp::numRowsub(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nb[i];
}

int OcpQp::numColsub(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return 1;
}

int OcpQp::numRowsC(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->ng[i];
}

int OcpQp::numColsC(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nx[i];
}

int OcpQp::numRowsD(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->ng[i];
}

int OcpQp::numColsD(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->nu[i];
}

int OcpQp::numRowslg(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->ng[i];
}

int OcpQp::numColslg(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return 1;
}

int OcpQp::numRowsug(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return dimensions->ng[i];
}

int OcpQp::numColsug(int i) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    return 1;
}

void OcpQp::setQ(int i, double *Q, int num_rows, int num_cols) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_rows = dimensions->nx[i], expected_num_cols = dimensions->nx[i];
    if (expected_num_rows != num_rows || expected_num_cols != num_cols)
        throw std::invalid_argument("Matrix has wrong size (" + std::to_string(num_rows) + " x " +
            std::to_string(num_cols) + "), expected size (" + std::to_string(expected_num_rows) +
            " x " + std::to_string(expected_num_cols) + ").");
    
    setQ(i, Q);
}

void OcpQp::setQ(int i, double *Q) {
    int num_rows = dimensions->nx[i], num_cols = dimensions->nx[i];
    int row_offset = dimensions->nu[i], col_offset = dimensions->nu[i];
    blasfeo_pack_dmat(num_rows, num_cols, Q, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
}

void OcpQp::setS(int i, double *S, int num_rows, int num_cols) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    // In HPIPM, S is an nu times nx matrix
    int expected_num_rows = dimensions->nu[i], expected_num_cols = dimensions->nx[i];
    if (expected_num_rows != num_rows || expected_num_cols != num_cols)
        throw std::invalid_argument("Matrix has wrong size (" + std::to_string(num_rows) + " x " +
            std::to_string(num_cols) + "), expected size (" + std::to_string(expected_num_rows) +
            " x " + std::to_string(expected_num_cols) + ").");
    
    setS(i, S);
}

void OcpQp::setS(int i, double *S) {
    // In HPIPM, S is an nu times nx matrix
    int num_rows = dimensions->nu[i], num_cols = dimensions->nx[i];
    int row_offset = dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, S, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
}

void OcpQp::setR(int i, double *R, int num_rows, int num_cols) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_rows = dimensions->nu[i], expected_num_cols = dimensions->nu[i];
    if (expected_num_rows != num_rows || expected_num_cols != num_cols)
        throw std::invalid_argument("Matrix has wrong size (" + std::to_string(num_rows) + " x " +
            std::to_string(num_cols) + "), expected size (" + std::to_string(expected_num_rows) +
            " x " + std::to_string(expected_num_cols) + ").");
    
    setR(i, R);
}

void OcpQp::setR(int i, double *R) {
    int num_rows = dimensions->nu[i], num_cols = dimensions->nu[i];
    int row_offset = 0, col_offset = 0;
    blasfeo_pack_dmat(num_rows, num_cols, R, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
}

void OcpQp::setq(int i, double *q, int num_elems) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_elems = dimensions->nx[i];
    if (expected_num_elems != num_elems)
        throw std::invalid_argument("Vector has wrong size (" + std::to_string(num_elems) + " x " +
            "1), expected size (" + std::to_string(expected_num_elems) + " x 1).");
    
    setq(i, q);
}

void OcpQp::setq(int i, double *q) {
    int num_rows = dimensions->nx[i], num_cols = 1;
    int row_offset = dimensions->nu[i] + dimensions->nx[i], col_offset = dimensions->nu[i];
    blasfeo_pack_tran_dmat(num_rows, num_cols, q, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
    blasfeo_pack_dvec(num_rows, q, &(qp->rq[i]), col_offset);
}

void OcpQp::setr(int i, double *r, int num_elems) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_elems = dimensions->nu[i];
    if (expected_num_elems != num_elems)
        throw std::invalid_argument("Vector has wrong size (" + std::to_string(num_elems) + " x " +
            "1), expected size (" + std::to_string(expected_num_elems) + " x 1).");
    
    setr(i, r);
}

void OcpQp::setr(int i, double *r) {
    int num_rows = dimensions->nu[i], num_cols = 1;
    int row_offset = dimensions->nu[i] + dimensions->nx[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, r, num_rows, &(qp->RSQrq[i]), row_offset, col_offset);
    blasfeo_pack_dvec(num_rows, r, &(qp->rq[i]), 0);
}

void OcpQp::setA(int i, double *A, int num_rows, int num_cols) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    int expected_num_rows = dimensions->nx[i+1], expected_num_cols = dimensions->nx[i];
    if (expected_num_rows != num_rows || expected_num_cols != num_cols)
        throw std::invalid_argument("Matrix has wrong size (" + std::to_string(num_rows) + " x " +
            std::to_string(num_cols) + "), expected size (" + std::to_string(expected_num_rows) +
            " x " + std::to_string(expected_num_cols) + ").");
    
    setA(i, A);
}

void OcpQp::setA(int i, double *A) {
    int num_rows = dimensions->nx[i+1], num_cols = dimensions->nx[i];
    int row_offset = dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, A, num_rows, &(qp->BAbt[i]), row_offset, col_offset);
}

void OcpQp::setB(int i, double *B, int num_rows, int num_cols) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    int expected_num_rows = dimensions->nx[i+1], expected_num_cols = dimensions->nu[i];
    if (expected_num_rows != num_rows || expected_num_cols != num_cols)
        throw std::invalid_argument("Matrix has wrong size (" + std::to_string(num_rows) + " x " +
            std::to_string(num_cols) + "), expected size (" + std::to_string(expected_num_rows) +
            " x " + std::to_string(expected_num_cols) + ").");
    
    setB(i, B);
}

void OcpQp::setB(int i, double *B) {
    int num_rows = dimensions->nx[i+1], num_cols = dimensions->nu[i];
    int row_offset = 0, col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, B, num_rows, &(qp->BAbt[i]), row_offset, col_offset);
}

void OcpQp::setb(int i, double *b, int num_elems) {
    if (i < 0 || i >= N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "[.");
    int expected_num_elems = dimensions->nx[i+1];
    if (expected_num_elems != num_elems)
        throw std::invalid_argument("Vector has wrong size (" + std::to_string(num_elems) + " x " +
            "1), expected size (" + std::to_string(expected_num_elems) + " x 1).");
    
    setb(i, b);
}

void OcpQp::setb(int i, double *b) {
    int num_rows = dimensions->nx[i+1], num_cols = 1;
    int row_offset = dimensions->nx[i] + dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, b, num_rows, &(qp->BAbt[i]), row_offset, col_offset);
    blasfeo_pack_dvec(num_rows, b, &(qp->b[i]), 0);
}

void OcpQp::setlb(int i, double *lb, int num_elems) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_elems = dimensions->nb[i];
    if (expected_num_elems != num_elems)
        throw std::invalid_argument("Vector has wrong size (" + std::to_string(num_elems) + " x " +
            "1), expected size (" + std::to_string(expected_num_elems) + " x 1).");
    
    setlb(i, lb);
}

void OcpQp::setlb(int i, double *lb) {
    int num_elems = dimensions->nb[i];
    int offset = 0;
    blasfeo_pack_dvec(num_elems, lb, &(qp->d[i]), offset);
}

void OcpQp::setub(int i, double *ub, int num_elems) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_elems = dimensions->nb[i];
    if (expected_num_elems != num_elems)
        throw std::invalid_argument("Vector has wrong size (" + std::to_string(num_elems) + " x " +
            "1), expected size (" + std::to_string(expected_num_elems) + " x 1).");
    
    setub(i, ub);
}

void OcpQp::setub(int i, double *ub) {
    int num_elems = dimensions->nb[i];
    int offset = dimensions->nb[i] + dimensions->ng[i];
    blasfeo_pack_dvec(num_elems, ub, &(qp->d[i]), offset);
    blasfeo_dvecsc(num_elems, -1.0, &(qp->d[i]), offset);
}

void OcpQp::setC(int i, double *C, int num_rows, int num_cols) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_rows = dimensions->nx[i], expected_num_cols = dimensions->ng[i];
    if (expected_num_rows != num_rows || expected_num_cols != num_cols)
        throw std::invalid_argument("Matrix has wrong size (" + std::to_string(num_rows) + " x " +
            std::to_string(num_cols) + "), expected size (" + std::to_string(expected_num_rows) +
            " x " + std::to_string(expected_num_cols) + ").");
    
    setC(i, C);
}

void OcpQp::setC(int i, double *C) {
    int num_rows = dimensions->nx[i], num_cols = dimensions->ng[i];
    int row_offset = dimensions->nu[i], col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, C, num_rows, &(qp->DCt[i]), row_offset, col_offset);
}

void OcpQp::setD(int i, double *D, int num_rows, int num_cols) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_rows = dimensions->nu[i], expected_num_cols = dimensions->ng[i];
    if (expected_num_rows != num_rows || expected_num_cols != num_cols)
        throw std::invalid_argument("Matrix has wrong size (" + std::to_string(num_rows) + " x " +
            std::to_string(num_cols) + "), expected size (" + std::to_string(expected_num_rows) +
            " x " + std::to_string(expected_num_cols) + ").");
    
    setD(i, D);
}

void OcpQp::setD(int i, double *D) {
    int num_rows = dimensions->nu[i], num_cols = dimensions->ng[i];
    int row_offset = 0, col_offset = 0;
    blasfeo_pack_tran_dmat(num_rows, num_cols, D, num_rows, &(qp->DCt[i]), row_offset, col_offset);
}

void OcpQp::setlg(int i, double *lg, int num_elems) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_elems = dimensions->ng[i];
    if (expected_num_elems != num_elems)
        throw std::invalid_argument("Vector has wrong size (" + std::to_string(num_elems) + " x " +
            "1), expected size (" + std::to_string(expected_num_elems) + " x 1).");
    
    setlg(i, lg);
}

void OcpQp::setlg(int i, double *lg) {
    int num_elems = dimensions->ng[i];
    int offset = dimensions->nb[i];
    blasfeo_pack_dvec(num_elems, lg, &(qp->d[i]), offset);
}

void OcpQp::setug(int i, double *ug, int num_elems) {
    if (i < 0 || i > N)
        throw std::out_of_range("Index " + std::to_string(i) + " needs to be in [0, N=" +
            std::to_string(N) + "].");
    int expected_num_elems = dimensions->ng[i];
    if (expected_num_elems != num_elems)
        throw std::invalid_argument("Vector has wrong size (" + std::to_string(num_elems) + " x " +
            "1), expected size (" + std::to_string(expected_num_elems) + " x 1).");
    
    setug(i, ug);
}

void OcpQp::setug(int i, double *ug) {
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


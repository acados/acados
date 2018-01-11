#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>

#include "acados_cpp/ocp_qp.h"

#include "acados_c/ocp_qp.h"
#include "acados/utils/print.h"

std::ostream& operator<<(std::ostream& oss, const OcpQp& qp) {
    print_ocp_qp_in(qp.qp);
    return oss;
}

void OcpQp::copyDimensions(std::vector<int> dimensions, int *dimension_ptr) {
    if (dimensions.size() != this->dimensions->N+1)
        throw std::invalid_argument("Dimensions must be defined for all N+1 nodes.");
    if (std::any_of(dimensions.begin(), dimensions.end(), [](int a) {return a < 0;}))
        throw std::invalid_argument("Dimension must be non-negative for all N+1 nodes.");
    std::copy(dimensions.begin(), dimensions.end(), dimension_ptr);
}

OcpQp::OcpQp(int N, std::vector<int> nx, std::vector<int> nu, std::vector<int> nb,
             std::vector<int> nc) :
    dimensions(nullptr),
    qp(nullptr)
{
    if (N <= 0) throw std::invalid_argument("Number of stages must be positive");

    dimensions = create_ocp_qp_dims(N);
    dimensions->N = N;

    copyDimensions(nx, dimensions->nx);
    copyDimensions(nu, dimensions->nu);
    copyDimensions(nb, dimensions->nb);
    copyDimensions(nc, dimensions->ng);

    qp = create_ocp_qp_in(dimensions);
}

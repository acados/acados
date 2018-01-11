#include <algorithm>
#include <iterator>

#include "acados_cpp/ocp_qp.h"

#include "acados_c/ocp_qp.h"


void OcpQp::setStateDimensionAt(int index, int numStates) {
    if (numStates < 0)
        throw std::invalid_argument("Number of states must be non-negative.");
    nx.at(index) = numStates;
    dimensions->nx[index] = numStates;
}

void OcpQp::setStateDimensions(int numStates) {
    if (numStates < 0)
        throw std::invalid_argument("Number of states must be non-negative.");
    std::fill(nx.begin(), nx.end(), numStates);
    std::fill(dimensions->nx, dimensions->nx + N+1, numStates);
}

void OcpQp::setControlDimensionAt(int index, int numControls) {
    if (numControls < 0)
        throw std::invalid_argument("Number of controls must be non-negative.");
    nu.at(index) = numControls;
    dimensions->nu[index] = numControls;
}

void OcpQp::setControlDimensions(int numControls) {
    if (numControls < 0)
        throw std::invalid_argument("Number of controls must be non-negative.");
    std::fill(nu.begin(), nu.end(), numControls);
    std::fill(dimensions->nu, dimensions->nu + N+1, numControls);
}

OcpQp::OcpQp(int N, int numStates, int numControls, int numBounds, int numConstraints) :
    N(N),
    nx(N),
    nu(N),
    nb(N),
    nc(N),
    dimensions(nullptr),
    qp(nullptr)
{
    dimensions = create_ocp_qp_dims(N);

    setStateDimensions(numStates);
    setControlDimensions(numControls);
    // By default, no controls on last stage
    setControlDimensionAt(0, 0);
    // If not specified, nx bounds on first stage
    if (numBounds == 0)
        setStateDimensionAt(0, numStates);

    qp = create_ocp_qp_in(dimensions);

    std::fill(nb.begin(), nb.end(), 0);
    std::fill(nc.begin(), nc.end(), 0);
}

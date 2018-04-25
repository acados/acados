
#include "acados_cpp/ocp_qp/ocp_qp_dimensions.hpp"

#include <algorithm>

using std::map;
using std::string;
using std::vector;

namespace acados {

/*
 * OCP QP dimensions should have fields
 * nx  - array with the number of states, per stage
 * nu  - array ""       ""        controls
 * 
 * and can have following optional fields
 * nbx - array ""       ""        state bounds
 * nbu - array ""       ""        control bounds
 * ng  - array ""       ""        polytopic constraints
 */
bool are_valid_ocp_qp_dimensions(const map<string, vector<uint>> dimensions) {

    if (!dimensions.count("nx") || dimensions.at("nx").size() == 0)
        return false;
    if (!dimensions.count("nu") || dimensions.at("nu").size() == 0)
        return false;

    if (dimensions.count("nu") && dimensions.at("nu").back() != 0)
        return false;
    if (dimensions.count("nbu") && dimensions.at("nbu").back() != 0)
        return false;

    std::vector<string> dim_names {"nx", "nu", "nbx", "nbu", "ng", "ns"};

    int expected_size = dimensions.at("nx").size();
    for (auto elem : dimensions) {
        if (std::find(dim_names.begin(), dim_names.end(), elem.first) == dim_names.end())
            return false;  // dimension name not valid
        if (elem.second.size() != expected_size)
            return false;
    }

    for (int i = 0; i < expected_size; ++i) {
        if (dimensions.at("nbx").at(i) > dimensions.at("nx").at(i))
            return false;
        if (dimensions.at("nbu").at(i) > dimensions.at("nu").at(i))
            return false;
    }
    return true;
}

std::unique_ptr<ocp_qp_dims> create_ocp_qp_dimensions_ptr(const map<string, vector<uint>>& dims) {

    if (!are_valid_ocp_qp_dimensions(dims))
        throw std::invalid_argument("Invalid dimensions map.");

    int N = dims.at("nx").size() - 1;

    auto dimensions_ptr = std::unique_ptr<ocp_qp_dims>(ocp_qp_dims_create(N));

    // required fields
    std::copy_n(std::begin(dims.at("nx")), N+1, dimensions_ptr->nx);
    std::copy_n(std::begin(dims.at("nu")), N+1, dimensions_ptr->nu);
    
    // optional fields
    vector<uint> nbx = dims.count("nbx") ? dims.at("nbx") : vector<uint>(N+1, 0);
    vector<uint> nbu = dims.count("nbu") ? dims.at("nbu") : vector<uint>(N+1, 0);
    vector<uint> ng = dims.count("ng") ? dims.at("ng") : vector<uint>(N+1, 0);
    vector<uint> ns = dims.count("ns") ? dims.at("ns") : vector<uint>(N+1, 0);

    std::vector<uint> nb(N+1);
    std::transform(nbx.begin(), nbx.end(), nbu.begin(), nb.begin(), std::plus<uint>());    

    std::copy_n(std::begin(nbx), N+1, dimensions_ptr->nbx);
    std::copy_n(std::begin(nbu), N+1, dimensions_ptr->nbu);
    std::copy_n(std::begin(ng), N+1, dimensions_ptr->ng);
    std::copy_n(std::begin(nb), N+1, dimensions_ptr->nb);

    return dimensions_ptr;

}

}  // namespace acados

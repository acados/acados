
#include "acados_cpp/ocp_dimensions.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

using std::map;
using std::string;
using std::vector;

namespace acados
{
bool are_valid_ocp_dimensions(const map<string, vector<int>>& dimensions,
                              const vector<string> valid_names)
{
    if (!dimensions.count("nx") || dimensions.at("nx").size() == 0) return false;
    if (!dimensions.count("nu") || dimensions.at("nu").size() == 0) return false;

    if (dimensions.count("nu") && dimensions.at("nu").back() != 0) return false;
    if (dimensions.count("nbu") && dimensions.at("nbu").back() != 0) return false;

    auto expected_size = dimensions.at("nx").size();
    for (auto elem : dimensions)
    {
        if (std::find(valid_names.begin(), valid_names.end(), elem.first) == valid_names.end())
            return false;  // dimension name not valid
        if (elem.second.size() != expected_size) return false;
    }

    for (size_t i = 0; i < expected_size; ++i)
    {
        if (dimensions.at("nbx").at(i) > dimensions.at("nx").at(i)) return false;
        if (dimensions.at("nbu").at(i) > dimensions.at("nu").at(i)) return false;
    }
    return true;
}

std::unique_ptr<ocp_qp_dims> create_ocp_qp_dimensions_ptr(const map<string, vector<int>>& dims)
{
    if (!are_valid_ocp_dimensions(dims, {"nx", "nu", "nbx", "nbu", "nb", "ng", "ns"}))
        throw std::invalid_argument("Invalid dimensions map.");

    int N = dims.at("nx").size() - 1;

    auto dimensions_ptr = std::unique_ptr<ocp_qp_dims>(ocp_qp_dims_create(N));

    // required fields
    std::copy_n(std::begin(dims.at("nx")), N + 1, dimensions_ptr->nx);
    std::copy_n(std::begin(dims.at("nu")), N + 1, dimensions_ptr->nu);

    // optional fields
    vector<int> nbx = dims.count("nbx") ? dims.at("nbx") : vector<int>(N + 1, 0);
    vector<int> nbu = dims.count("nbu") ? dims.at("nbu") : vector<int>(N + 1, 0);
    vector<int> ng = dims.count("ng") ? dims.at("ng") : vector<int>(N + 1, 0);
    vector<int> ns = dims.count("ns") ? dims.at("ns") : vector<int>(N + 1, 0);

    std::vector<int> nb(N + 1);
    std::transform(nbx.begin(), nbx.end(), nbu.begin(), nb.begin(), std::plus<int>());

    std::copy_n(std::begin(nbx), N + 1, dimensions_ptr->nbx);
    std::copy_n(std::begin(nbu), N + 1, dimensions_ptr->nbu);
    std::copy_n(std::begin(ng), N + 1, dimensions_ptr->ng);
    std::copy_n(std::begin(nb), N + 1, dimensions_ptr->nb);

    return dimensions_ptr;
}

// std::unique_ptr<ocp_nlp_dims>
// create_ocp_nlp_dimensions_ptr(const map<string, vector<int>>& dims) {

//     if (!are_valid_ocp_dimensions(dims, {"nx", "nu", "nbx", "nbu", "nb", "ng", "nh", "ns"}))
//         throw std::invalid_argument("Invalid dimensions map.");

//     int N = dims.at("nx").size() - 1;

//     auto dimensions_ptr = std::unique_ptr<ocp_nlp_dims>(ocp_nlp_dims_create(N));

//     // required fields
//     std::copy_n(std::begin(dims.at("nx")), N+1, dimensions_ptr->nx);
//     std::copy_n(std::begin(dims.at("nu")), N+1, dimensions_ptr->nu);

//     // optional fields
//     vector<int> nbx = dims.count("nbx") ? dims.at("nbx") : vector<int>(N+1, 0);
//     vector<int> nbu = dims.count("nbu") ? dims.at("nbu") : vector<int>(N+1, 0);
//     vector<int> ng = dims.count("ng") ? dims.at("ng") : vector<int>(N+1, 0);
//     vector<int> nh = dims.count("nh") ? dims.at("nh") : vector<int>(N+1, 0);
//     vector<int> ns = dims.count("ns") ? dims.at("ns") : vector<int>(N+1, 0);

//     // calculate nv
//     std::transform(nx.begin(), nx.end(), nu.begin(), nv.begin(), std::plus<int>());
//     std::transform(nv.begin(), nv.end(), ns.begin(), nv.begin(), std::plus<int>());

//     // calculate ni
//     std::transform(nbx.begin(), nbx.end(), nbu.begin(), ni.begin(), std::plus<int>());
//     std::transform(ni.begin(), ni.end(), ng.begin(), ni.begin(), std::plus<int>());
//     std::transform(ni.begin(), ni.end(), nh.begin(), ni.begin(), std::plus<int>());

//     std::copy_n(std::begin(nv), N+1, dimensions_ptr->nv);
//     std::copy_n(std::begin(ni), N+1, dimensions_ptr->ni);

//     return dimensions_ptr;

// }

}  // namespace acados

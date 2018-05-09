
#include "acados_cpp/ocp_bounds.hpp"

#include <algorithm>
#include <stdexcept>

namespace acados
{
using std::vector;

vector<int> calculate_idxb(const vector<double>& lower_bound, const vector<double>& upper_bound)
{
    if (lower_bound.size() != upper_bound.size())
        throw std::invalid_argument("Lower bound must have same shape as upper bound.");

    vector<int> bound_indices;

    for (size_t idx = 0; idx < lower_bound.size(); ++idx)
        if (lower_bound.at(idx) != -INFINITY || upper_bound.at(idx) != +INFINITY)
            bound_indices.push_back(idx);  // there is a one-sided or two-sided bound at this index

    return bound_indices;
}

vector<vector<int>> calculate_all_idxb(const vector<vector<double>>& lower_bounds,
                                       const vector<vector<double>>& upper_bounds)
{
    if (lower_bounds.size() != upper_bounds.size())
        throw std::invalid_argument("Number of lower and upper bounds must be the same.");

    vector<vector<int>> bounds_indices(lower_bounds.size());

    std::transform(lower_bounds.begin(), lower_bounds.end(), upper_bounds.begin(),
                   bounds_indices.begin(),
                   [](vector<double> lb, vector<double> ub) { return calculate_idxb(lb, ub); });

    return bounds_indices;
}

void copy_at(vector<double>& output, vector<double> input, vector<int> copy_ids)
{
    if (copy_ids.size() > output.size())
        std::invalid_argument("Given more copy indices than expected length of output");

    for (auto idx : copy_ids) output.at(idx) = input.at(idx);
}

}  // namespace acados

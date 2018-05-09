
#ifndef INTERFACES_ACADOS_CPP_OCP_BOUNDS_HPP_
#define INTERFACES_ACADOS_CPP_OCP_BOUNDS_HPP_

#include <cmath>
#include <vector>

namespace acados
{
/*
 * From lower and upper bounds, calculate the bounded indices.
 */
std::vector<int> calculate_idxb(const std::vector<double>& lower_bound,
                                const std::vector<double>& upper_bound);

/*
 * From a collection of lower bounds and upper bounds, return all bounded indices.
 */
std::vector<std::vector<int>> calculate_all_idxb(
    const std::vector<std::vector<double>>& lower_bounds,
    const std::vector<std::vector<double>>& upper_bounds);

/*
 * Copy values from input at given indices, pad remaining elements with a given value.
 */
void copy_at(std::vector<double>& output, std::vector<double> input, std::vector<int> copy_ids);

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_BOUNDS_HPP_

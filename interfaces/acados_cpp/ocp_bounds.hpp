
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OCP_BOUNDS_HPP_H
#define ACADOS_INTERFACES_ACADOS_CPP_OCP_BOUNDS_HPP_H

#include <cmath>
#include <vector>

/*
 * From lower and upper bounds, calculate the bounded indices.
 */
std::vector<int>
calculate_idxb(const std::vector<double>& lower_bound, const std::vector<double>& upper_bound);

/*
 * From a collection of lower bounds and upper bounds, return all bounded indices.
 */
std::vector<std::vector<int>>
calculate_all_idxb(const std::vector<std::vector<double>>& lower_bounds, const std::vector<std::vector<double>>& upper_bounds);

/*
 * 
 */
void
copy_at(std::vector<double>& output, std::vector<double> input, std::vector<int> copy_ids);

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OCP_BOUNDS_HPP_H

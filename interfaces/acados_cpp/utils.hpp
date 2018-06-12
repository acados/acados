
#ifndef INTERFACES_ACADOS_CPP_UTILS_HPP_
#define INTERFACES_ACADOS_CPP_UTILS_HPP_

#include <string>
#include <utility>
#include <vector>

namespace acados
{
std::string to_string(std::pair<int, int> p);

template <typename T>
std::string to_string(std::vector<T> v);

bool match(std::pair<int, int> dims, int nb_elems);

template <typename T>
const T& clamp(const T& lo, const T& hi, const T& val);

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_UTILS_HPP_

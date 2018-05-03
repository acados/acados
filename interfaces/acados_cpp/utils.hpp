
#ifndef ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_

#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <utility>
#include <string>
#include <vector>

namespace acados {

std::string to_string(std::pair<int, int> p);

template<typename T>
std::string to_string(std::vector<T> v);

bool match(std::pair<int, int> dims, int nb_elems);

template<typename T>
const T& clamp(const T& lo, const T& hi, const T& val);

std::string load_error_message();

}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_
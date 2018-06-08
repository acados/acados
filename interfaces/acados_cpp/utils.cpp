
#include "acados_cpp/utils.hpp"

namespace acados
{
std::string to_string(std::pair<int, int> p)
{
    return "( " + std::to_string(p.first) + ", " + std::to_string(p.second) + " )";
}

template <typename T>
std::string to_string(std::vector<T> v)
{
    std::string result_string = " vector of length " + std::to_string(v.size()) + ": [\n ";
    for (auto it : v)
    {
        result_string += std::to_string(it) + ", ";
    }
    return result_string + "]\n";
}

bool match(std::pair<int, int> dims, int nb_elems)
{
    int nb_expected_elems = dims.first * dims.second;
    if (nb_expected_elems == 0 || nb_expected_elems == nb_elems) return true;
    return false;
}

template <typename T>
const T& clamp(const T& lo, const T& hi, const T& val)
{
    if (val < lo)
        return lo;
    else if (val > hi)
        return hi;

    return val;
}

}  // namespace acados

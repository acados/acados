
#ifndef ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_

#include <utility>
#include <vector>

namespace std {
    std::string to_string(std::pair<uint, uint> p) {
        return "( " + std::to_string(p.first) + ", " + std::to_string(p.second) + " )";
    }

    template<typename T>
    std::string to_string(std::vector<T> v) {
        std::string result_string = " vector of length " + std::to_string(v.size()) + ": [\n ";
        for(auto it : v) {
            result_string += std::to_string(it) + ", ";
        }
        return result_string + "]\n";
    }
}  // namespace std

#endif  // ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_
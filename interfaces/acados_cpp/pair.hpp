
#ifndef ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_

namespace std {
    std::string to_string(std::pair<uint, uint> p) {
        return "( " + std::to_string(p.first) + ", " + std::to_string(p.second) + " )";
    }
}  // namespace std

#endif  // ACADOS_INTERFACES_ACADOS_CPP_PAIR_HPP_
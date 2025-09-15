#ifndef {{ package_name | upper }}_UTILS_HPP
#define {{ package_name | upper }}_UTILS_HPP

#include <ostream>
#include <array>
#include <algorithm>

template <size_t N>
std::ostream& operator<<(std::ostream& os, const double (&arr)[N]) {
    os << "[";
    for (size_t i = 0; i < N; ++i) {
        os << arr[i];
        if (i + 1 < N) os << ", ";
    }
    os << "]";
    return os;
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << "[";
    for (size_t i = 0; i < N; ++i) {
        os << arr[i];
        if (i < N - 1) os << ", ";
    }
    os << "]";
    return os;
}

#endif // {{ package_name | upper }}_UTILS_HPP

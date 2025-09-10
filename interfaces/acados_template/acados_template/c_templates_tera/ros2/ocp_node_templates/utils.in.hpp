#ifndef {{ ros_opts.package_name | upper }}_UTILS_HPP
#define {{ ros_opts.package_name | upper }}_UTILS_HPP

#include <rclcpp/rclcpp.hpp>
#include <ostream>
#include <vector>
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

/**
 * @brief Extract the diagonal of a square matrix stored as a flat row-major std::array.
 *
 * @tparam T element type
 * @tparam N matrix dimension (NxN)
 * @param flat_mat flat row-major storage of size N*N
 * @return std::array<T, N> diagonal elements
 */
template<typename T, size_t N>
inline std::array<T, N> diagonal(const std::array<T, N * N>& flat_mat) noexcept
{
    std::array<T, N> d{};
    for (size_t i = 0; i < N; ++i) {
        d[i] = flat_mat[i * N + i];
    }
    return d;
}

/**
 * @brief Create a diagonal matrix (flat row-major) from a vector.
 *
 * @tparam T element type
 * @tparam N vector dimension (N)
 * @param v flat vector storage of size N
 * @return std::array<T, N * N> with only diagonal elements
 */
template<typename T, size_t N>
inline std::array<T, N * N> diag_from_vec(const std::array<T, N>& v) noexcept
{
    std::array<T, N * N> mat{};
    // mat is zero-initialized; set only diagonal entries
    for (size_t i = 0; i < N; ++i) {
        mat[i * N + i] = v[i];
    }
    return mat;
}

/**
 * @brief Create an alternating array, which takes a pattern of two values and expands it to a larger size.
 *
 * @tparam T element type
 * @tparam N vector dimension (N)
 * @param pattern_values flat pattern vector storage of size 2
 * @return std::array<T, N>
 */
template<typename T, size_t N, size_t K>
inline std::array<T, N> create_repeating_array(const std::array<T, K>& pattern_values) {
    std::array<T, N> alternated{};
    for (size_t i = 0; i < N; ++i) {
        alternated[i] = pattern_values[i % K];
    }
    return alternated;
}


inline std::vector<int> range(int start, int end) {
    std::vector<int> result;
    if (end <= start) return result;
    result.reserve(end - start);
    for (int i = start; i < end; ++i) {
        result.push_back(i);
    }
    return result;
}

#endif // {{ ros_opts.package_name | upper }}_UTILS_HPP

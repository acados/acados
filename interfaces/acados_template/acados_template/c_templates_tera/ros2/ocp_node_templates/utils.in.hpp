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

template <size_t N>
inline void get_and_check_array_param(rclcpp::Node* node, const std::string& param_name, std::array<double, N>& destination) {
    auto param_value = node->get_parameter(param_name).as_double_array();

    if (param_value.size() != N) {
        RCLCPP_ERROR(node->get_logger(), "Parameter '%s' has the wrong size. Expected: %ld, got: %ld",
                     param_name.c_str(), N, param_value.size());
        return;
    }
    std::copy_n(param_value.begin(), N, destination.begin());
}

template <size_t N>
inline void update_param_array(const rclcpp::Parameter& param, 
                        std::array<double, N>& destination_array,
                        rcl_interfaces::msg::SetParametersResult& result)
{
    auto values = param.as_double_array();

    if (values.size() != N) {
        result.successful = false;
        result.reason = "Parameter '" + param.get_name() + "' has size " +
                        std::to_string(values.size()) + ", but expected is " + std::to_string(N) + ".";
        return;
    }

    std::copy_n(values.begin(), N, destination_array.begin());
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

#endif // {{ ros_opts.package_name | upper }}_UTILS_HPP
#ifndef TEST_TEST_UTILS_EIGEN_H_
#define TEST_TEST_UTILS_EIGEN_H_

/* Ignore compiler warnings from Eigen */
#if defined(__clang__)
    #include <Eigen/Dense>
#elif defined(__GNUC__)
    #if __GNUC__ > 4
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wmisleading-indentation"
        #include <Eigen/Dense>
        #pragma GCC diagnostic pop
    #else
        #include <Eigen/Dense>
    #endif
#endif

#endif  // TEST_TEST_UTILS_EIGEN_H_

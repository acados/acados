#ifndef TEST_TEST_UTILS_EIGEN_H_
#define TEST_TEST_UTILS_EIGEN_H_

/* Ignore compiler warnings from Eigen */
#if defined(__clang__)
    #include <Eigen/Dense>
#elif defined(__GNUC__)
    #if __GNUC__ == 6
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wmisleading-indentation"
        #include <Eigen/Dense>
        #pragma GCC diagnostic pop
    #endif
#endif

#endif  // TEST_TEST_UTILS_EIGEN_H_

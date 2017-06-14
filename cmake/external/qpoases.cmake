include(ExternalProject)

ExternalProject_Add(
    qpoases_project

    CONFIGURE_COMMAND ""
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/qpOASES"
    BINARY_DIR "${PROJECT_SOURCE_DIR}/external/qpOASES/build"
    CONFIGURE_COMMAND cmake -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ..
    BUILD_COMMAND cmake --build .
    # INSTALL_COMMAND cmake --build . --target install
    INSTALL_COMMAND ""
    LOG_CONFIGURE 1  # suppress output
    LOG_BUILD 1
)

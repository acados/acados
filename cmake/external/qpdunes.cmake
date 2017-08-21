include(ExternalProject)

ExternalProject_Add(
    qpdunes_project

    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/qpDUNES-dev"
    BINARY_DIR "${PROJECT_SOURCE_DIR}/external/qpDUNES-dev/build"
    CONFIGURE_COMMAND cmake ..
    BUILD_COMMAND cmake --build .
    INSTALL_COMMAND ""
    LOG_CONFIGURE 1  # suppress output
    LOG_BUILD 1
)

ExternalProject_Get_Property(qpdunes_project source_dir)

add_library(qpdunes STATIC IMPORTED GLOBAL)
add_dependencies(qpdunes qpdunes_project)
set_property(TARGET qpdunes PROPERTY IMPORTED_LOCATION "${source_dir}/build/lib/libqpdunes.a")

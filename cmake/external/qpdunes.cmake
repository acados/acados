include(ExternalProject)

ExternalProject_Add(
    qpdunes_project

    CONFIGURE_COMMAND ""
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/qpDUNES-dev"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(qpdunes_project source_dir)

add_library(qpdunes STATIC IMPORTED)
add_dependencies(qpdunes qpdunes_project)
set_property(TARGET qpdunes PROPERTY IMPORTED_LOCATION "${source_dir}/build/lib/libqpdunes.a")

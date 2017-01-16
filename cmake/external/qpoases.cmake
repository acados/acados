include(ExternalProject)

ExternalProject_Add(
    qpoases_project

    CONFIGURE_COMMAND ""
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/qpoases"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(qpoases_project source_dir)
add_library(qpoases STATIC IMPORTED)
add_dependencies(qpoases qpoases_project)
set_property(TARGET qpoases PROPERTY IMPORTED_LOCATION "${source_dir}/bin/libqpOASES_e.a")

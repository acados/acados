include(ExternalProject)

ExternalProject_Add(
    ma27_project

    CONFIGURE_COMMAND ./configure
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/coinhsl"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    # TODO(dimitris): maybe not with the install command
    INSTALL_COMMAND cp ${PROJECT_SOURCE_DIR}/external/coinhsl/.libs/libcoinhsl.a ${PROJECT_SOURCE_DIR}/external/OOQP/libma27.a
)

ExternalProject_Get_Property(ma27_project source_dir)

# add_library(ma27 STATIC IMPORTED)
# add_dependencies(ma27 ma27_project)
# set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION "${project_source_dir}/external/OOQP/libma27.a")

include(ExternalProject)

ExternalProject_Add(
    ma27_project

    CONFIGURE_COMMAND ./configure
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/coinhsl"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)

ExternalProject_Add_Step(ma27_project cp_lib
   COMMAND cp ${PROJECT_SOURCE_DIR}/external/coinhsl/.libs/libcoinhsl.a ${PROJECT_SOURCE_DIR}/external/OOQP/libma27.a
   DEPENDEES build
)

ExternalProject_Get_Property(ma27_project source_dir)
add_library(ma27 STATIC IMPORTED)
add_dependencies(ma27 ma27_project)
set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/OOQP/libma27.a")

include(ExternalProject)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Needed for cross-compiling
    set(HOST_FLAG "--host=${COMPILER_PREFIX}")
endif()

ExternalProject_Add(
    ma27_project

    CONFIGURE_COMMAND ./configure "${HOST_FLAG}" "CC=${CMAKE_C_COMPILER}" "CFLAGS=-O2 -fPIC" "FCFLAGS=-O2 -fPIC"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/coinhsl"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Add_Step(ma27_project cp_lib
   COMMAND cp ${PROJECT_SOURCE_DIR}/external/coinhsl/.libs/libcoinhsl.a ${PROJECT_SOURCE_DIR}/external/OOQP/libma27.a
   DEPENDEES build
)

ExternalProject_Get_Property(ma27_project source_dir)
add_library(ma27 STATIC IMPORTED GLOBAL)
add_dependencies(ma27 ma27_project)
set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/OOQP/libma27.a")

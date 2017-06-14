include(ExternalProject)

find_package(OpenBLAS REQUIRED)
find_package(FortranLibs REQUIRED)
include(external/ma27)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(OOQP_LDFLAGS "-lc++")
else()
    set(OOQP_LDFLAGS "")
endif()

ExternalProject_Add(
    ooqp_project

    DEPENDS ma27
    CONFIGURE_COMMAND ./configure "CXXFLAGS=-O2 -fPIC" "CFLAGS=-O2 -fPIC" "FFLAGS=-O2 -fPIC" "LDFLAGS=${OOQP_LDFLAGS}"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/OOQP"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND ""
    LOG_CONFIGURE 1  # suppress output
    LOG_BUILD 1
)

ExternalProject_Get_Property(ooqp_project source_dir)

add_library(ooqpgensparse STATIC IMPORTED GLOBAL)
add_dependencies(ooqpgensparse ooqp_project)
set_property(TARGET ooqpgensparse PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpgensparse.a")

add_library(ooqpsparse STATIC IMPORTED GLOBAL)
add_dependencies(ooqpsparse ooqp_project)
set_property(TARGET ooqpsparse PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpsparse.a")

add_library(ooqpgondzio STATIC IMPORTED GLOBAL)
add_dependencies(ooqpgondzio ooqp_project)
set_property(TARGET ooqpgondzio PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpgondzio.a")

add_library(ooqpbase STATIC IMPORTED GLOBAL)
add_dependencies(ooqpbase ma27 ooqp_project)
set_property(TARGET ooqpbase PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpbase.a")

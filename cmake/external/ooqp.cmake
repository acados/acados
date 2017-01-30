include(ExternalProject)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(OOQP_LDFLAGS "-lc++")
else()
    set(OOQP_LDFLAGS " ")
endif()

ExternalProject_Add(
    ooqp_project

    CONFIGURE_COMMAND ./configure
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/OOQP"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make LDFLAGS="${OOQP_LDFLAGS}"
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(ooqp_project source_dir)

add_dependencies(ooqp_project ma27_project)

add_library(ooqpgensparse STATIC IMPORTED)
add_dependencies(ooqpgensparse ooqp_project)
set_property(TARGET ooqpgensparse PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpgensparse.a")

add_library(ooqpsparse STATIC IMPORTED)
add_dependencies(ooqpsparse ooqp_project)
set_property(TARGET ooqpsparse PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpsparse.a")

add_library(ooqpgondzio STATIC IMPORTED)
add_dependencies(ooqpgondzio ooqp_project)
set_property(TARGET ooqpgondzio PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpgondzio.a")

add_library(ooqpbase STATIC IMPORTED)
add_dependencies(ooqpbase ma27 ooqp_project)
set_property(TARGET ooqpbase PROPERTY IMPORTED_LOCATION "${source_dir}/lib/libooqpbase.a")

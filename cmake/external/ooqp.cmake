include(ExternalProject)

find_package(OpenBLAS REQUIRED)
add_library(openblas UNKNOWN IMPORTED)
set_property(TARGET openblas PROPERTY IMPORTED_LOCATION ${OpenBLAS_LIB})

find_package(FortranLibs REQUIRED)
add_library(gfortran UNKNOWN IMPORTED)
set_property(TARGET gfortran PROPERTY IMPORTED_LOCATION ${FORTRAN_LIBRARY})

include(external/ma27)

set(OOQP_LDFLAGS "")
set(HOST_FLAG "")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(OOQP_LDFLAGS "-lc++")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Needed for cross-compiling
    set(HOST_FLAG "--host=${COMPILER_PREFIX}")
endif()

ExternalProject_Add(
    ooqp_project

    DEPENDS ma27
    PREFIX "${PROJECT_BINARY_DIR}/external/ooqp"
    CONFIGURE_COMMAND ./configure "--prefix=${PROJECT_BINARY_DIR}/external/ooqp/" "${HOST_FLAG}"
                                  "CXX=${CMAKE_CXX_COMPILER}" "CXXFLAGS=-O2 -fPIC" "CC=${CMAKE_C_COMPILER}"
                                  "CFLAGS=-O2 -fPIC" "FFLAGS=-O2 -fPIC" "LDFLAGS=${OOQP_LDFLAGS}"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/OOQP"
    BUILD_COMMAND make clean all
    INSTALL_COMMAND make install
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Get_Property(ooqp_project BINARY_DIR)

ExternalProject_Add_Step(ooqp_project copy_ooqp
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${EXTERNAL_SRC_DIR}/ooqp/" "${BINARY_DIR}"
    DEPENDERS configure
)

ExternalProject_Add_Step(ooqp_project create_lib_folder
    COMMAND ${CMAKE_COMMAND} -E chdir ${BINARY_DIR} mkdir lib
    DEPENDERS configure
)

ExternalProject_Add_Step(ooqp_project copy_ma27
    COMMAND cp ${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a ${BINARY_DIR}/libma27.a
    DEPENDERS configure
)

add_library(ooqp INTERFACE)
target_link_libraries(ooqp INTERFACE
    ooqpgensparse
    ooqpsparse
    ooqpgondzio
    ooqpbase
    ma27
    openblas
    gfortran
    c++
    m)

set_property(TARGET ooqp
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/external/ooqp/include>
        $<INSTALL_INTERFACE:include>)

install(DIRECTORY ${BINARY_DIR}/include/
    DESTINATION include/ooqp
    FILES_MATCHING PATTERN "*.h")

add_library(ooqpgensparse STATIC IMPORTED GLOBAL)
add_dependencies(ooqpgensparse ooqp_project)
set_property(TARGET ooqpgensparse PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpgensparse.a")
install(FILES "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpgensparse.a" DESTINATION lib)

add_library(ooqpsparse STATIC IMPORTED GLOBAL)
add_dependencies(ooqpsparse ooqp_project)
set_property(TARGET ooqpsparse PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpsparse.a")
install(FILES "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpsparse.a" DESTINATION lib)

add_library(ooqpgondzio STATIC IMPORTED GLOBAL)
add_dependencies(ooqpgondzio ooqp_project)
set_property(TARGET ooqpgondzio PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpgondzio.a")
install(FILES "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpgondzio.a" DESTINATION lib)

add_library(ooqpbase STATIC IMPORTED GLOBAL)
add_dependencies(ooqpbase ma27 ooqp_project)
set_property(TARGET ooqpbase PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpbase.a")
install(FILES "${PROJECT_BINARY_DIR}/external/ooqp/lib/libooqpbase.a" DESTINATION lib)

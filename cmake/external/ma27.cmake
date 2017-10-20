include(ExternalProject)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Needed for cross-compiling
    set(HOST_FLAG "--host=${COMPILER_PREFIX}")
endif()

ExternalProject_Add(
    ma27_project
    PREFIX "${PROJECT_BINARY_DIR}/external/ma27/"
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND sh "${EXTERNAL_SRC_DIR}/coinhsl/configure" "--prefix=${PROJECT_BINARY_DIR}/external/ma27/" "${HOST_FLAG}" "CC=${CMAKE_C_COMPILER}" "CFLAGS=-O2 -fPIC" "FCFLAGS=-O2 -fPIC"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/coinhsl"
    BUILD_COMMAND make clean all
    INSTALL_COMMAND make install
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Add_Step(ma27_project rename_library
    COMMAND mv ${PROJECT_BINARY_DIR}/external/ma27/lib/libcoinhsl.a ${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a
    DEPENDEES install
)

add_library(ma27 STATIC IMPORTED GLOBAL)
add_dependencies(ma27 ma27_project)
set_property(TARGET ma27 PROPERTY IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a")
install(FILES "${PROJECT_BINARY_DIR}/external/ma27/lib/libma27.a" DESTINATION lib)

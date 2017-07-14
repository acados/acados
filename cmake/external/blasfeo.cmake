include(ExternalProject)

if(NOT DEFINED BLASFEO_TARGET)
#   set(BLASFEO_TARGET X64_INTEL_HASWELL)
    # set(BLASFEO_TARGET X64_INTEL_SANDY_BRIDGE)
  set(BLASFEO_TARGET GENERIC)
endif()
if(NOT DEFINED BLASFEO_LA)
#    set(BLASFEO_LA HIGH_PERFORMANCE)
   set(BLASFEO_LA REFERENCE)
#   set(BLASFEO_LA BLAS)
endif()

set(OS "")
set(RANLIB_COMMAND "")
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(OS "LINUX")
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(OS "MAC")
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(OS "WINDOWS")
    set(RANLIB_COMMAND ${CMAKE_RANLIB} libblasfeo.a)
endif()

ExternalProject_Add(
    blasfeo_project

    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/blasfeo"
    BINARY_DIR "${PROJECT_SOURCE_DIR}/external/blasfeo/build"
    CONFIGURE_COMMAND cmake -DOS=${OS} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DTARGET=${BLASFEO_TARGET} -DLA=${BLASFEO_LA} -DCMAKE_VERBOSE_MAKEFILE=ON ..
    BUILD_COMMAND cmake --build . --target blasfeo -- OS=${OS} CC=${CMAKE_C_COMPILER} TARGET=${BLASFEO_TARGET} LA=${BLASFEO_LA}
    INSTALL_COMMAND "${RANLIB_COMMAND}"
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Get_Property(blasfeo_project source_dir)
add_library(blasfeo STATIC IMPORTED GLOBAL)
add_dependencies(blasfeo blasfeo_project)
set_property(TARGET blasfeo PROPERTY IMPORTED_LOCATION "${source_dir}/build/libblasfeo.a")

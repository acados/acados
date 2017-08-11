include(ExternalProject)

if(NOT DEFINED HPMPC_TARGET)
#	set(HPMPC_TARGET X64_AVX2)
#   set(HPMPC_TARGET X64_AVX)
	set(HPMPC_TARGET C99_4X4)
endif()
if(NOT DEFINED HPMPC_BLASFEO_PATH)
	set(HPMPC_BLASFEO_PATH ${PROJECT_SOURCE_DIR}/external/blasfeo)
endif()

set(OS "")
set(RANLIB_COMMAND "")
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(OS "LINUX")
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(OS "MAC")
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(OS "WINDOWS")
    set(RANLIB_COMMAND ${CMAKE_RANLIB} libhpmpc.a)
endif()

ExternalProject_Add(
    hpmpc_project

    DEPENDS blasfeo
    CONFIGURE_COMMAND make clean
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/hpmpc"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_MAKE_COMMAND} static_library -j 2 OS=${OS} CC=${CMAKE_C_COMPILER} TARGET=${HPMPC_TARGET} BLASFEO_PATH=${HPMPC_BLASFEO_PATH}
    INSTALL_COMMAND "${RANLIB_COMMAND}"
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Get_Property(hpmpc_project source_dir)
add_library(hpmpc STATIC IMPORTED GLOBAL)
add_dependencies(hpmpc hpmpc_project blasfeo_project)
set_property(TARGET hpmpc PROPERTY IMPORTED_LOCATION "${source_dir}/libhpmpc.a")

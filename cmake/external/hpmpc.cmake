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
    set(RANLIB_COMMAND ${CMAKE_RANLIB} ${PROJECT_BINARY_DIR}/external/hpmpc/lib/libhpmpc.a)  # for mingw-w64
endif()

ExternalProject_Add(
    hpmpc_project

    DEPENDS blasfeo
    PREFIX "${PROJECT_BINARY_DIR}/external/hpmpc"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/hpmpc"
    CONFIGURE_COMMAND ls && make clean
    BUILD_COMMAND make static_library -j 2 OS=${OS} CC=${CMAKE_C_COMPILER} TARGET=${HPMPC_TARGET} BLASFEO_PATH=${HPMPC_BLASFEO_PATH}
    INSTALL_COMMAND make install_static PREFIX=${PROJECT_BINARY_DIR}/external
    # LOG_CONFIGURE 1  # suppress output
    # LOG_BUILD 1
)

ExternalProject_Get_Property(hpmpc_project BINARY_DIR)

ExternalProject_Add_Step(hpmpc_project copy_hpmpc
    COMMAND ${CMAKE_COMMAND} -E copy_directory "${EXTERNAL_SRC_DIR}/hpmpc/" "${BINARY_DIR}"
    DEPENDERS configure
)

ExternalProject_Add_Step(hpmpc_project ranlib_step
    COMMAND ${RANLIB_COMMAND}
    DEPENDEES install
)

ExternalProject_Add_Step(hpmpc_project copy_target_h
    COMMAND cp ${BINARY_DIR}/include/target.h ${EXTERNAL_SRC_DIR}/hpmpc/include/target.h
    DEPENDEES install
)

add_library(hpmpc STATIC IMPORTED GLOBAL)
add_dependencies(hpmpc hpmpc_project blasfeo_project)
set_property(TARGET hpmpc PROPERTY IMPORTED_LOCATION "${BINARY_DIR}/libhpmpc.a")

install(FILES "${BINARY_DIR}/libhpmpc.a" DESTINATION lib)

install(DIRECTORY ${EXTERNAL_SRC_DIR}/hpmpc/include/
    DESTINATION include/hpmpc
    FILES_MATCHING PATTERN "*.h")

set_property(TARGET hpmpc
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/external/>
        $<INSTALL_INTERFACE:include/hpmpc>)

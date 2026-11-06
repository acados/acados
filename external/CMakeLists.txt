# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error -w")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error -w")
endif()

if(CMAKE_SYSTEM_NAME MATCHES "dSpace")
    set(USE_C99_MATH OFF CACHE BOOL "fmax is not implemented for dSpace.")
endif()

# When BLASFEO is used as part of acados
if(ACADOS_WITH_SYSTEM_BLASFEO)
    find_package(blasfeo REQUIRED)
else()
    set(BLASFEO_HEADERS_INSTALLATION_DIRECTORY "include/blasfeo/include" CACHE STRING "")

    set(TARGET ${BLASFEO_TARGET} CACHE STRING "Set CPU architecture target" FORCE)
    set(MF "PANELMAJ" CACHE STRING "Matrix format" FORCE)
    set(BLAS_API OFF CACHE BOOL "Compile BLAS API" FORCE)
    add_subdirectory(blasfeo)
endif()

# When HPIPM is used as part of acados
if(ACADOS_WITH_SYSTEM_BLASFEO)
    set(HPIPM_FIND_BLASFEO ON CACHE BOOL "" FORCE)
else()
    set(HPIPM_BLASFEO_PATH ${BLASFEO_SRC_DIR} CACHE STRING "Set CPU architecture target" FORCE)
endif()
set(HPIPM_BLASFEO_LIB OFF CACHE BOOL "Turn off BLASFEO from HPIPM" FORCE)
set(HPIPM_HEADERS_INSTALLATION_DIRECTORY "include/hpipm/include" CACHE STRING "")
set(TARGET ${HPIPM_TARGET} CACHE STRING "Set CPU architecture target" FORCE)
set(HPIPM_TESTING OFF CACHE BOOL "Tests disabled" FORCE)
add_subdirectory(hpipm)

if(ACADOS_WITH_HPMPC)
	# TODO: The following definition should go ASAP!
	# When hpmpc is used as part of acados
	set(HPMPC_HEADERS_INSTALLATION_DIRECTORY "include/hpmpc/include" CACHE STRING "")
	add_subdirectory(hpmpc)
endif()

if(ACADOS_WITH_QORE)
    set(QORE_DENSE TRUE CACHE BOOL "Build dense qore library")
    add_subdirectory(qore)
endif()

if(ACADOS_WITH_QPDUNES)
    add_subdirectory(qpdunes)
endif()

if(ACADOS_WITH_QPOASES)
    add_subdirectory(qpoases)
endif()

if(ACADOS_WITH_DAQP)
    set(PROFILING OFF CACHE BOOL "Turn off profiling in DAQP")
    set(SOFT_WEIGHTS ON CACHE BOOL "Turn on soft weights in DAQP")
    add_subdirectory(daqp)
endif()

if(ACADOS_WITH_OSQP)
    set(OSQP_ENABLE_PROFILING OFF CACHE BOOL "Turn off profiling in OSQP")
    set(OSQP_ENABLE_INTERRUPT OFF CACHE BOOL "Turn off user interrupt in OSQP")
    add_subdirectory(osqp)
endif()


if(ACADOS_WITH_CLARABEL)
    set(CLARABEL_BUILD_EXAMPLES OFF CACHE BOOL "set Clarabel examples" FORCE)
    # The Clarabel library is built by cargo and copied to this directory in the build tree,
    # which exists independently of the install prefix.
    set(CLARABEL_C_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/Clarabel.cpp" CACHE STRING "Set Clarabel.cpp output dir" FORCE)
    set(PROFILING OFF CACHE BOOL "Turn off profiling in CLARABEL")
    set(CTRLC OFF CACHE BOOL "Turn off CTRLC in CLARABEL")
    add_subdirectory(Clarabel.cpp)
    add_dependencies(acados libclarabel_c_shared)

    # libclarabel_c_shared is part of the exported link interface of acados.
    # Its include directories (source tree) and library path (build tree) are only valid
    # when building acados. For the installed export, link against the library installed
    # to <prefix>/lib. No installed acados header requires clarabel.h.
    get_target_property(CLARABEL_INCLUDE_DIRS libclarabel_c_shared INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(CLARABEL_LINK_LIBS libclarabel_c_shared INTERFACE_LINK_LIBRARIES)
    set(CLARABEL_INCLUDE_DIRS_EXPORT "")
    foreach(dir IN LISTS CLARABEL_INCLUDE_DIRS)
        list(APPEND CLARABEL_INCLUDE_DIRS_EXPORT "$<BUILD_INTERFACE:${dir}>")
    endforeach()
    set(CLARABEL_LINK_LIBS_EXPORT "")
    foreach(lib IN LISTS CLARABEL_LINK_LIBS)
        string(FIND "${lib}" "${CLARABEL_C_OUTPUT_DIR}/" CLARABEL_LIB_POS)
        if(CLARABEL_LIB_POS EQUAL 0)
            get_filename_component(CLARABEL_LIB_NAME "${lib}" NAME)
            list(APPEND CLARABEL_LINK_LIBS_EXPORT
                "$<BUILD_INTERFACE:${lib}>"
                "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/lib/${CLARABEL_LIB_NAME}>")
        else()
            list(APPEND CLARABEL_LINK_LIBS_EXPORT "${lib}")
        endif()
    endforeach()
    set_target_properties(libclarabel_c_shared PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CLARABEL_INCLUDE_DIRS_EXPORT}"
        INTERFACE_LINK_LIBRARIES "${CLARABEL_LINK_LIBS_EXPORT}")
endif()

if (ACADOS_WITH_OOQP)
    include(external/ooqp)
endif()

# Build list with all source files to go into the acados library
file(GLOB_RECURSE ACADOS_SRC ${PROJECT_SOURCE_DIR}/*.c)
# Exclude helper files
list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/condensing_helper_functions.c)
# Exclude main.c file for compiling on ESP32
list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/main.c)

if (NOT EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/ocp_qp_ooqp.c)
endif()

# Sources for examples
#set(HPMPC_EXAMPLE_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_hpmpc.c)
#set(CONDENSING_QPOASES_EXAMPLE_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_condensing_qpoases.c)
set(CHEN_MODEL_SRC ${PROJECT_SOURCE_DIR}/Chen_model.c)
set(NMPC_EXAMPLE_SRC ${PROJECT_SOURCE_DIR}/test_nmpc.c)

# Merge options from HPMPC Cmake
if(NOT DEFINED HPMPC_TARGET)
#	set(HPMPC_TARGET X64_AVX2)
	set(HPMPC_TARGET X64_AVX)
#	set(HPMPC_TARGET C99_4X4)
endif()

# Merge options from BLASFEO Cmake
if(NOT DEFINED BLASFEO_TARGET)
#   set(BLASFEO_TARGET X64_INTEL_HASWELL)
    set(BLASFEO_TARGET X64_INTEL_SANDY_BRIDGE)
#   set(BLASFEO_TARGET GENERIC)
endif()
if(NOT DEFINED BLASFEO_LA)
    set(BLASFEO_LA HIGH_PERFORMANCE)
#   set(BLASFEO_LA REFERENCE)
#   set(BLASFEO_LA BLAS)
endif()

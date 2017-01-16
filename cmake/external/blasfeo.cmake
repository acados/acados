include(ExternalProject)

if(NOT DEFINED BLASFEO_TARGET)
#	set(BLASFEO_TARGET X64_INTEL_HASWELL)
	set(BLASFEO_TARGET X64_INTEL_SANDY_BRIDGE)
#	set(BLASFEO_TARGET GENERIC)
endif()
if(NOT DEFINED BLASFEO_LA)
	set(BLASFEO_LA HIGH_PERFORMANCE)
#   set(BLASFEO_LA REFERENCE)
#	set(BLASFEO_LA BLAS)
endif()

ExternalProject_Add(
    blasfeo_project

    CONFIGURE_COMMAND ""
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/blasfeo"
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make clean static_library -j 2 TARGET=${BLASFEO_TARGET} LA=${BLASFEO_LA}
    INSTALL_COMMAND ""
)

ExternalProject_Get_Property(blasfeo_project source_dir)
add_library(blasfeo STATIC IMPORTED)
add_dependencies(blasfeo blasfeo_project)
set_property(TARGET blasfeo PROPERTY IMPORTED_LOCATION "${source_dir}/libblasfeo.a")

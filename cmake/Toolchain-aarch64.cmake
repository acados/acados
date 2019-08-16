set(CMAKE_SYSTEM_NAME "Generic" CACHE STRING "")
message(CMAKE_C_COMPILER)

#### user settings start ###
set(COMPILER_PATH "/opt/Xilinx/SDK/2018.3/gnu/aarch64/lin/aarch64-none/bin")
set(CROSS_PREFIX   "${COMPILER_PATH}/aarch64-none-elf-" CACHE STRING "")

link_directories("path-to-bsp/lib")
include_directories("path-to-bsp/include")
#e.g. include_directories("/home/username/project/a53_1_bsp/psu_cortexa53_1/include")

set(LINKER_SCRIPT_PATH 
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform/lscript_ultrazed.ld" 
	CACHE STRING "Linker script file")
#### user settings end ###

set(HPIPM_TARGET "GENERIC" CACHE STRING "HPIPM Target architecture")
set(BLASFEO_TARGET "ARMV8A_ARM_CORTEX_A53" CACHE STRING "BLASFEO Target architecture")

set (CMAKE_SYSTEM_PROCESSOR "aarch64"            CACHE STRING "")

set(CMAKE_C_FLAGS_DEBUG "-O0 -g3" CACHE STRING "")
set(CMAKE_C_FLAGS_RELEASE "-O2" CACHE STRING "")

set (CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER CACHE STRING "")
set (CMAKE_FIND_ROOT_PATH_MODE_LIBRARY NEVER CACHE STRING "")
set (CMAKE_FIND_ROOT_PATH_MODE_INCLUDE NEVER CACHE STRING "")

set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set (CMAKE_C_COMPILER   "${CROSS_PREFIX}gcc" CACHE STRING "")
set (CMAKE_CXX_COMPILER "${CROSS_PREFIX}g++" CACHE STRING "")

set(EMBEDDED_TARGET "XILINX_NONE_ELF" CACHE STRING "Xilinx bare-metal")
add_definitions(-D__XILINX_NONE_ELF__)
remove_definitions(-DLINUX)
remove_definitions(-D__LINUX__)

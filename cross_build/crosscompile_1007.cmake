# !! ANDREA: RUN THIS FORM A WINDOWS SHELL TO AVOID PATH FORMAT PROBLEMS. THEN COMPILE FROM CYGWIN. !!
# 1) run cmake -DCMAKE_TOOLCHAIN_FILE="crosscompile_1007.cmake" -G "Unix Makefiles" -DACADOS_WITH_OSQP=OFF 
#       -DACADOS_WITH_QPDUNES=OFF -DACADOS_WITH_HPMPC=OFF -DCMAKE_INSTALL_PREFIX="..\acados_install_dir" ..
# 
# 2) set the Windows environment variables QNX_TARGET and QNX_HOST flags to QNX_TARGET = 
# "C:\Program Files\dSPACE RCPHIL 2017-A\Compiler\QNX650_520\target\qnx6" and
#       QNX_HOST = 
#       "C:\Program Files\dSPACE RCPHIL 2017-A\Compiler\QNX650_520\host\win32\x86" 
# - when using cmake with Git Bash: -DCMAKE_SH="CMAKE_SH-NOTFOUND"

set(CMAKE_MAKE_PROGRAM "C:/Program\ Files/dSPACE RCPHIL\ 2017-A/Exe/DSMAKE.exe")
set(CMAKE_C_COMPILER   "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-g++-5.2.0.exe") #  # ntoppc-gcc-5.2.0.exe
set(CMAKE_CXX_COMPILER "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-g++-5.2.0.exe") #  # ntoppc-gcc-5.2.0.exe
set(CMAKE_ASM_COMPILER "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-as.exe") #  # ntoppc-gcc-5.2.0.exe
set(CMAKE_ASM_COMPILER "C:/Program\ Files/dSPACE\ RCPHIL\ 2017-A/Compiler/QNX650_520/host/win32/x86/usr/bin/ntoppc-ar.exe") #  # ntoppc-gcc-5.2.0.exe

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x c" CACHE STRING "" FORCE) # needed when using ntoppcg++ as done in dspace makefile
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mpowerpc -mcpu=e500mc  -mcpu=e500mc -mtune=e500mc64 -mhard-float -mdouble-float -EB " CACHE STRING "" FORCE)

set(CMAKE_C_COMPILER_WORKS 1 CACHE INTERNAL "")
set(CMAKE_CXX_COMPILER_WORKS 1 CACHE INTERNAL "")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


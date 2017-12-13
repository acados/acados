# the name of the target operating system
set(CMAKE_SYSTEM_NAME Windows)

# Choose an appropriate compiler prefix
# for 32 or 64 bits mingw-w64
# see http://mingw-w64.sourceforge.net/
set(BITNESS "x86_64")
set(COMPILER_PREFIX "${BITNESS}-w64-mingw32")

# which compilers to use for C and C++
find_program(CMAKE_RC_COMPILER NAMES ${COMPILER_PREFIX}-windres)
find_program(CMAKE_C_COMPILER NAMES ${COMPILER_PREFIX}-gcc)
find_program(CMAKE_CXX_COMPILER NAMES ${COMPILER_PREFIX}-g++)
find_program(CMAKE_AR NAMES ${COMILER-PREFIX}-ar)
find_program(CMAKE_RANLIB NAMES ${COMILER-PREFIX}-ranlib)

# here is the target environment located
set(CMAKE_FIND_ROOT_PATH
    "/usr/local/opt/mingw-w64/toolchain-${BITNESS}/${COMPILER_PREFIX}/"
    "/usr/lib/gcc/${COMPILER_PREFIX}/*"
    "$ENV{HOME}/WindowsLibs"
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

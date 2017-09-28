set(CMAKE_SYSTEM_NAME dSpace)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

file(TO_CMAKE_PATH "C:\\Program Files\\dSPACE RCPHIL 2016-B" DSPACE_TOOLS)
set(DSPACE_PPCTOOLS ${DSPACE_TOOLS}/Compiler/PPCTools)
find_program(CMAKE_C_COMPILER NAMES ${DSPACE_PPCTOOLS}/bin/mccppc.exe)
find_program(CMAKE_CXX_COMPILER NAMES ${DSPACE_PPCTOOLS}/bin/cccppc.exe)
find_program(CMAKE_ASM_COMPILER NAMES ${DSPACE_PPCTOOLS}/bin/asmppc.exe)
find_program(CMAKE_AR NAMES ${DSPACE_PPCTOOLS}/bin/cccppc.exe)

set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> <OBJECTS> -o <TARGET>")
set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> <OBJECTS> -o <TARGET>")

find_program(CMAKE_MAKE_PROGRAM NAMES ${DSPACE_TOOLS}/Exe/DSMAKE.EXE)

set(CMAKE_FIND_ROOT_PATH ${DSPACE_PPCTOOLS})

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

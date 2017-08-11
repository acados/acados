message(STATUS "DSPACE")

set(CMAKE_STATIC_LIBRARY_PREFIX "" CACHE STRING "Static library prefix")
set(CMAKE_STATIC_LIBRARY_SUFFIX ".lib" CACHE STRING "Static library suffix")

set(CMAKE_INCLUDE_FLAG_C "-J")
set(CMAKE_INCLUDE_FLAG_CXX "-J")

add_definitions(-D__DSPACE__)
remove_definitions(-DLINUX)
remove_definitions(-D__LINUX__)

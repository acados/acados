
# The dSpace compiler only works with the prefixes/suffixes below!

set(CMAKE_IMPORT_LIBRARY_PREFIX "")
set(CMAKE_SHARED_LIBRARY_PREFIX "")
set(CMAKE_SHARED_MODULE_PREFIX  "")
set(CMAKE_STATIC_LIBRARY_PREFIX "")

set(CMAKE_EXECUTABLE_SUFFIX     ".exe")
set(CMAKE_IMPORT_LIBRARY_SUFFIX ".lib")
set(CMAKE_SHARED_LIBRARY_SUFFIX ".lib")
set(CMAKE_SHARED_MODULE_SUFFIX  ".lib")
set(CMAKE_STATIC_LIBRARY_SUFFIX ".lib")

set(CMAKE_C_FLAGS "\"-J${DSPACE_RTLIB}\"")

set(CMAKE_INCLUDE_FLAG_C "-J")
set(CMAKE_INCLUDE_FLAG_CXX "-J")

add_definitions(-D__DSPACE__)
remove_definitions(-DLINUX)
remove_definitions(-D__LINUX__)

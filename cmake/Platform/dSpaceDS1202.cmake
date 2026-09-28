#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.



# The dSpace compiler only works with the prefixes/suffixes below!

SET(CMAKE_SHARED_LIBRARY_PREFIX "lib")
SET(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
SET(CMAKE_STATIC_LIBRARY_PREFIX "lib")
SET(CMAKE_STATIC_LIBRARY_SUFFIX ".a")
SET(BUILD_SHARED_LIBS "OFF")

add_definitions(-DWINDOWS_SKIP_PTR_ALIGNMENT_CHECK)
remove_definitions(-DLINUX)
remove_definitions(-D__LINUX__)

#
#    This file is part of acados.
#
#    acados is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 3 of the License, or (at your option) any later version.
#
#    acados is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with acados; if not, write to the Free Software Foundation,
#    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#

include(external/blasfeo)
include(external/hpmpc)
add_subdirectory(qpDUNES-dev)
add_subdirectory(qpOASES)
if (NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/OOQP)
        include(external/ooqp)
    endif()
endif()

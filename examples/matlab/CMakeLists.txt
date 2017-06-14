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

find_package(Matlab)
if(SWIG_MATLAB AND MATLAB_FOUND)
add_test(NAME ocp_qp_matlab
    COMMAND ${MATLAB_EXECUTABLE} -noFigureWindows -nodesktop -nosplash -r
        "try, ocp_qp_example, catch ME, disp(ME.getReport()), exit(1), end, exit(0)")

add_test(NAME ocp_nlp_matlab
    COMMAND ${MATLAB_EXECUTABLE} -noFigureWindows -nodesktop -nosplash -r
        "try, ocp_nlp_example, catch ME, disp(ME.getReport()), exit(1), end, exit(0)")
endif()

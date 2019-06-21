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

find_package(Matlab COMPONENTS MAIN_PROGRAM REQUIRED)

add_test(NAME ocp_qp_matlab
    # run() is not useful here as it changes cwd, so we have to use addpath()
    COMMAND ${Matlab_MAIN_PROGRAM} "-noFigureWindows" "-nodesktop" "-nosplash" "-nosoftwareopengl" "-r"
        "addpath('${CMAKE_CURRENT_SOURCE_DIR}'); try, ocp_qp_example, catch ME, disp(ME.getReport()), exit(1), end, exit(0)")

# add_test(NAME ocp_nlp_matlab
#     COMMAND ${Matlab_MAIN_PROGRAM} "-noFigureWindows" "-nodesktop" "-nosplash" "-nosoftwareopengl" "-r"
#         "addpath('${CMAKE_CURRENT_SOURCE_DIR}'); try, ocp_nlp_example, catch ME, disp(ME.getReport()), exit(1), end, exit(0)")

# add_test(NAME sim_matlab
    # COMMAND ${Matlab_MAIN_PROGRAM} "-noFigureWindows" "-nodesktop" "-nosplash" "-nosoftwareopengl" "-r"
        # "addpath('${CMAKE_CURRENT_SOURCE_DIR}'); try, sim_example, catch ME, disp(ME.getReport()), exit(1), end, exit(0)")

set_tests_properties(ocp_qp_matlab
    PROPERTIES
        ENVIRONMENT "MATLABPATH=${CMAKE_BINARY_DIR}/swig/matlab/:$ENV{MATLABPATH}")

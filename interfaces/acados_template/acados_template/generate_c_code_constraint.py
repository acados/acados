#   This file is part of acados.
#
#   acados is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 3 of the License, or (at your option) any later version.
#
#   acados is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with acados; if not, write to the Free Software Foundation,
#   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#

from casadi import *
import os

def generate_c_code_constraint( constraint, suffix_name ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if  casadi_version != '3.4.0':
        # old casadi versions
        raise Exception('Please download and install Casadi 3.4.0 to ensure compatibility with acados. Version ' + casadi_version + ' currently in use.')

    # load constraint variables and expression
    x = constraint.x
    u = constraint.u
    # nc = nh or np 
    nc = constraint.nc 
    con_exp = constraint.expr
    con_name = constraint.name

    # get dimensions
    nx = x.size()[0]
    nu = u.size()[0]

    # set up functions to be exported
    fun_name = con_name + suffix_name
    # TODO(andrea): first output seems to be ignored in the C code
    jac_x = jacobian(con_exp, x);
    jac_u = jacobian(con_exp, u);
    constraint_fun_jac_tran = Function(fun_name, [vertcat(u,x)], [con_exp, vertcat(transpose(jac_u), transpose(jac_x))])

    # generate C code
    if not os.path.exists('c_generated_code'):
        os.mkdir('c_generated_code')

    os.chdir('c_generated_code')
    gen_dir = con_name + suffix_name 
    if not os.path.exists(gen_dir):
        os.mkdir(gen_dir)
    gen_dir_location = './' + gen_dir
    os.chdir(gen_dir_location)
    file_name = con_name + suffix_name
    constraint_fun_jac_tran.generate(file_name, casadi_opts)
    os.chdir('../..')

    return

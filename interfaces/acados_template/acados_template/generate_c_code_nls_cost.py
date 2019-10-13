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

def function generate_c_code_nls_cost( cost, suffix_name ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if casadi_version not in ('3.4.5', '3.4.0'):
        # old casadi versions
        raise Exception('Please download and install CasADi 3.4.0 to ensure compatibility with acados. Version ' + casadi_version + ' currently in use.')

    # load cost variables and expression
    x = constraint.x
    u = constraint.u
    cost_exp = cost.expr
    cost_name = cost.name

    # get dimensions
    nx = x.size()[0]
    nu = u.size()[0]

    # set up functions to be exported
    fun_name = con_name + suffix_name

    cost_jac_exp = jacobian(cost_exp, vertcat([u, x]))
    
    nls_cost_fun = Function( fun_name, [x, u], \
            [ cost_exp, cost_jac_exp ])

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

    nls_cost_fun.generate( file_name, casadi_opts )
    
    os.chdir('../..')
    return


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

def generate_c_code_explicit_ode( model ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if  casadi_version != '3.4.0':
        # old casadi versions
        raise Exception('Please download and install Casadi 3.4.0 to ensure compatibility with acados. Version ' + casadi_version + ' currently in use.')

    # load model
    x = model.x
    u = model.u
    f_expl = model.f_expl_expr
    model_name = model.name

    ## get model dimensions
    nx = x.size()[0]
    nu = u.size()[0]

    ## set up functions to be exported
    Sx = SX.sym('Sx', nx, nx)
    Sp = SX.sym('Sp', nx, nu)
    lambdaX = SX.sym('lambdaX', nx, 1)

    fun_name = model_name + '_expl_ode_fun'
    expl_ode_fun = Function(fun_name, [x,u], [f_expl])
    # TODO: Polish: get rid of SX.zeros
    vdeX = SX.zeros(nx,nx)
    vdeX = vdeX + jtimes(f_expl,x,Sx)

    vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u)
    vdeP = vdeP + jtimes(f_expl,x,Sp)

    fun_name = model_name + '_expl_vde_forw'
    expl_vde_forw = Function(fun_name, [x,Sx,Sp,u], [f_expl,vdeX,vdeP])

    jacX = SX.zeros(nx,nx) + jacobian(f_expl,x)

    adj = jtimes(f_expl, vertcat(x, u), lambdaX, True)

    fun_name = model_name + '_expl_vde_adj'
    expl_vde_adj = Function(fun_name, [x,lambdaX,u], [adj])

    S_forw = vertcat(horzcat(Sx, Sp), horzcat(DM.zeros(nu,nx), DM.eye(nu)))
    hess = mtimes(transpose(S_forw),jtimes(adj, vertcat(x,u), S_forw))
    hess2 = []
    for j in range(nx+nu):
        for i in range(j,nx+nu):
            hess2 = vertcat(hess2, hess[i,j])

    fun_name = model_name + '_expl_ode_hess'
    expl_ode_hess = Function(fun_name, [x, Sx, Sp, lambdaX, u], [adj, hess2])

    ## generate C code
    if not os.path.exists('c_generated_code'):
        os.mkdir('c_generated_code')

    os.chdir('c_generated_code')
    model_dir = model_name + '_model'
    if not os.path.exists(model_dir):
        os.mkdir(model_dir)
    model_dir_location = './' + model_dir
    os.chdir(model_dir_location)
    fun_name = model_name + '_expl_ode_fun'
    expl_ode_fun.generate(fun_name, casadi_opts)

    fun_name = model_name + '_expl_vde_forw'
    expl_vde_forw.generate(fun_name, casadi_opts)
    
    fun_name = model_name + '_expl_vde_adj'
    expl_vde_adj.generate(fun_name, casadi_opts)
    
    fun_name = model_name + '_expl_ode_hess'
    expl_ode_hess.generate(fun_name, casadi_opts)
    os.chdir('../..')

    return

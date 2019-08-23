#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#


from casadi import *
import os

def generate_c_code_explicit_ode( model ):

    casadi_version = CasadiMeta.version()
    casadi_opts = dict(mex=False, casadi_int='int', casadi_real='double')

    if  casadi_version not in ('3.4.5', '3.4.0'):
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
    if str(type(f_expl)) == "<class 'casadi.casadi.SX'>":
        Sx = SX.sym('Sx', nx, nx)
        Sp = SX.sym('Sp', nx, nu)
        lambdaX = SX.sym('lambdaX', nx, 1)
    elif str(type(f_expl)) == "<class 'casadi.casadi.MX'>":
        Sx = MX.sym('Sx', nx, nx)
        Sp = MX.sym('Sp', nx, nu)
        lambdaX = MX.sym('lambdaX', nx, 1)
    else:
        raise Exception("Invalid type for f_expl! Possible types are 'SX' and 'MX'. Exiting.")

    fun_name = model_name + '_expl_ode_fun'
    expl_ode_fun = Function(fun_name, [x,u], [f_expl])
    # TODO: Polish: get rid of SX.zeros
    if str(type(f_expl)) == "<class 'casadi.casadi.SX'>":
        vdeX = SX.zeros(nx,nx)
    else: 
        vdeX = MX.zeros(nx,nx)

    vdeX = vdeX + jtimes(f_expl,x,Sx)

    if str(type(f_expl)) == "<class 'casadi.casadi.SX'>":
        vdeP = SX.zeros(nx,nu) + jacobian(f_expl,u)
    else: 
        vdeP = MX.zeros(nx,nu) + jacobian(f_expl,u)

    vdeP = vdeP + jtimes(f_expl,x,Sp)

    fun_name = model_name + '_expl_vde_forw'
    expl_vde_forw = Function(fun_name, [x,Sx,Sp,u], [f_expl,vdeX,vdeP])

    if str(type(f_expl)) == "<class 'casadi.casadi.SX'>":
        jacX = SX.zeros(nx,nx) + jacobian(f_expl,x)
    else: 
        jacX = MX.zeros(nx,nx) + jacobian(f_expl,x)

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

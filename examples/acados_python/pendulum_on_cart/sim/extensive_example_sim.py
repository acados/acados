#
# Copyright (c) The acados authors.
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

import sys, json
import time

sys.path.insert(0, '../common')

from acados_template import AcadosSim, AcadosSimSolver, AcadosModel, acados_dae_model_json_dump, sim_get_default_cmake_builder, GnsfModel
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np


import casadi as ca

def idx_perm_to_ipiv(idx_perm):
    n = len(idx_perm)
    vec = list(range(n))
    ipiv = np.zeros(n)

    for ii in range(n):
        idx0 = idx_perm[ii]
        for jj in range(ii,n):
            if vec[jj]==idx0:
                idx1 = jj
                break
        tmp = vec[ii]
        vec[ii] = vec[idx1]
        vec[idx1] = tmp
        ipiv[ii] = idx1

    ipiv = ipiv-1 # C 0-based indexing
    return ipiv


def export_pendulum_ode_model_with_gnsf_def(sim) -> AcadosModel:

    model_name = 'pendulum'

    # constants
    m_cart = 1. # mass of the cart [kg]
    m = 0.1 # mass of the ball [kg]
    g = 9.81 # gravity constant [m/s^2]
    l = 0.8 # length of the rod [m]

    # set up states & controls
    x1      = ca.MX.sym('x1')
    theta   = ca.MX.sym('theta')
    v1      = ca.MX.sym('v1')
    dtheta  = ca.MX.sym('dtheta')

    x = ca.vertcat(x1, theta, v1, dtheta)

    F = ca.MX.sym('F')
    u = ca.vertcat(F)

    # xdot
    x1_dot      = ca.MX.sym('x1_dot')
    theta_dot   = ca.MX.sym('theta_dot')
    v1_dot      = ca.MX.sym('v1_dot')
    dtheta_dot  = ca.MX.sym('dtheta_dot')

    xdot = ca.vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot)

    # parameters
    p = []

    # dynamics
    cos_theta = ca.cos(theta)
    sin_theta = ca.sin(theta)
    denominator = m_cart + m - m*cos_theta*cos_theta
    f_expl = ca.vertcat(v1,
                     dtheta,
                     (-m*l*sin_theta*dtheta*dtheta + m*g*cos_theta*sin_theta+F)/denominator,
                     (-m*l*cos_theta*sin_theta*dtheta*dtheta + F*cos_theta+(m_cart+m)*g*sin_theta)/(l*denominator)
                     )

    f_impl = xdot - f_expl
    model = AcadosModel()
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    # model.z = z
    model.p = p
    model.name = model_name

    # GNSF definition with x1 as linear output state
    E = np.eye(4)
    A = np.zeros((4,4))
    B = np.zeros((4,1))
    C = np.eye(4)
    phi = f_expl[:]
    y = ca.vertcat(x1, theta, v1, dtheta)
    uhat = u
    c = np.zeros((4,1))
    E_LO = np.eye(0)
    A_LO = np.zeros((0,0))
    B_LO = np.zeros((0,0))
    f_lo = 0
    c_LO = np.zeros((0,0))

    L_x = ca.jacobian(y, x)
    L_xdot = ca.jacobian(y, xdot)
    L_z = ca.MX.zeros(4,0)
    L_u = ca.jacobian(y, u)

    nontrivial_f_LO = 0
    purely_linear = 0
    idx_perm_x_1 = np.array([0, 1, 2, 3])
    ipiv_x = idx_perm_to_ipiv(idx_perm_x_1)
    print("ipiv_x:", ipiv_x)
    ipiv_z = np.array([])

    x1 = model.x[:]
    x1dot = model.xdot[:]
    z1 = ca.MX.sym("z1", 0, 0)

    model.gnsf_model = GnsfModel(
        x1=x1,
        x1dot=x1dot,
        z1=z1,
        y=y,
        uhat=uhat,
        phi=phi,
        f_LO=f_lo,
        E=E,
        A=A,
        B=B,
        C=C,
        c=c,
        E_LO=E_LO,
        A_LO=A_LO,
        B_LO=B_LO,
        c_LO=c_LO,
        L_x=L_x,
        L_xdot=L_xdot,
        L_z=L_z,
        L_u=L_u,
        nontrivial_f_LO=nontrivial_f_LO,
        purely_linear=purely_linear,
        # idx_perm_x_1=idx_perm_x_1,
        ipiv_x=ipiv_x,
        ipiv_z=ipiv_z,
    )

    # store meta information
    model.x_labels = ['$x$ [m]', r'$\theta$ [rad]', '$v$ [m]', r'$\dot{\theta}$ [rad/s]']
    model.u_labels = ['$F$ [N]']
    model.t_label = '$t$ [s]'

    return model


def main(gnsf_definition_mode):

    if gnsf_definition_mode not in ['user_provided', 'detected', 'imported']:
        raise ValueError("gnsf_definition_mode must be one of 'user_provided', 'detected', 'imported'")

    sim = AcadosSim()

    # model
    model = export_pendulum_ode_model()
    sim.model = model

    Tf = 0.1
    nx = model.x.rows()
    nu = model.u.rows()
    N = 200

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.num_stages = 7
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 10 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    sim.solver_options.integrator_type = "GNSF" # ERK, IRK, GNSF
    sim.solver_options.sens_forw = True
    sim.solver_options.sens_adj = True
    sim.solver_options.sens_hess = False
    sim.solver_options.sens_algebraic = False
    sim.solver_options.output_z = False
    sim.solver_options.sim_method_jac_reuse = False


    # if sim.solver_options.integrator_type == "GNSF":
    #     # Perform GNSF structure detection in Octave
    #     # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/external/casadi-octave
    #     # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/
    #     # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/acados_template_mex/

    #     # acados_dae_model_json_dump(model)
    #     # status = os.system('octave convert_dae2gnsf.m')
    #     with open(model.name + '_gnsf_functions.json', 'r') as f:
    #         gnsf_dict = json.load(f)
    #     sim.gnsf_model = gnsf_dict

    if sim.solver_options.integrator_type == "GNSF":
        if gnsf_definition_mode == 'imported':
            from acados_template import acados_dae_model_json_dump
            import os
            acados_dae_model_json_dump(model)
            # Set up Octave to be able to run the following:
            ## if using a virtual python env, the following lines can be added to the env/bin/activate script:
            # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/external/casadi-octave
            # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/
            # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/acados_template_mex/
            # echo
            # echo "OCTAVE_PATH=$OCTAVE_PATH"
            # status = os.system(
            #     "octave --eval \"convert_dae2gnsf({})\"".format("\'"+model.name+"_acados_dae.json\'")
            # )
            # if status == 0:
            #     print("\nsuccessfully detected GNSF structure in Octave\n")
            # else:
            #     Exception("Failed to detect GNSF structure in Octave")
            # load gnsf from json
            with open(model.name + '_gnsf_functions.json', 'r') as f:
                import json
                gnsf_dict = json.load(f)
            sim.gnsf_model = gnsf_dict
            breakpoint()
        elif gnsf_definition_mode == 'user_provided':
            sim.model = export_pendulum_ode_model_with_gnsf_def(sim)
        elif gnsf_definition_mode == 'detected':
            pass

    # create
    cmake_builder = sim_get_default_cmake_builder()
    acados_integrator = AcadosSimSolver(sim, cmake_builder=cmake_builder)

    simX = np.zeros((N+1, nx))
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    u0 = np.array([0.0])
    acados_integrator.set("u", u0)

    simX[0,:] = x0

    ## Single test call
    t0 = time.time()
    acados_integrator.set("seed_adj", np.ones((nx, 1)))
    acados_integrator.set("x", x0)
    acados_integrator.set("u", u0)
    status = acados_integrator.solve()
    time_external = time.time() - t0

    S_forw = acados_integrator.get("S_forw")
    Sx = acados_integrator.get("Sx")
    Su = acados_integrator.get("Su")
    S_adj = acados_integrator.get("S_adj")
    print(f"\ntimings of last call to acados_integrator: with Python interface, set and get {time_external*1e3:.4f}ms")


    # get timings (of last call)
    CPUtime = acados_integrator.get("CPUtime")
    LAtime = acados_integrator.get("LAtime")
    ADtime = acados_integrator.get("ADtime")
    print(f"\ntimings of last call to acados_integrator: overall CPU: {CPUtime*1e3:.4f} ms, linear algebra {LAtime*1e3:.4f} ms, external functions {ADtime*1e3:.4f} ms")

    print("S_forw, sensitivities of simulation result wrt x,u:\n", S_forw)
    print("Sx, sensitivities of simulation result wrt x:\n", Sx)
    print("Su, sensitivities of simulation result wrt u:\n", Su)
    print("S_adj, adjoint sensitivities:\n", S_adj)

    # turn off sensitivity propagation when not needed
    acados_integrator.options_set('sens_forw', False)
    acados_integrator.options_set('sens_adj', False)
    acados_integrator.options_set('sens_hess', False)

    # call in loop:
    for i in range(N):
        # set initial state
        acados_integrator.set("x", simX[i,:])
        # solve
        status = acados_integrator.solve()
        # get solution
        simX[i+1,:] = acados_integrator.get("x")

    if status != 0:
        raise Exception(f'acados returned status {status}.')


    # plot results
    plot_pendulum(np.linspace(0, N*Tf, N+1), 10, np.zeros((N, nu)), simX, latexify=False)


if __name__ == '__main__':
    # main(gnsf_definition_mode='imported')
    main(gnsf_definition_mode='user_provided')
    main(gnsf_definition_mode='detected')

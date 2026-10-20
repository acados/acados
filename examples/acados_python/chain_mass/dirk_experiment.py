#!/usr/bin/env python
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

# authors: Florian Messerer, Jonathan Frey

import os
import argparse
import numpy as np

from acados_template import AcadosSim, AcadosSimSolver
from export_chain_mass_model import export_chain_mass_model
from utils import *
from plot_utils import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("formulation", choices=["radau", "legendre", "esdirk", "sdirk"],
                        help="Which collocation type to use")
    parser.add_argument("num_stages", type=int,
                        help="number of collocationstages")
    parser.add_argument("--jac-reuse", action="store_true",
                        help="reuse jacobian")
    parser.add_argument("--exploit-dirk", action="store_true",
                        help="Exploit DIRK structure")
    parser.add_argument("--show-plots", action="store_true",
                        help="Show plots")
    parser.add_argument("-N", "--n-masses", type=int, default=10,
                        help="Number of masses")
    args = parser.parse_args()
    if args.formulation == "radau":
        COLLOC_TYPE = "GAUSS_RADAU_IIA"
    elif args.formulation == "legendre":
        COLLOC_TYPE = "GAUSS_LEGENDRE"
    elif args.formulation == "esdirk":
        COLLOC_TYPE = "ESDIRK"
    elif args.formulation == "sdirk":
        COLLOC_TYPE = "SDIRK"

    PLOT = args.show_plots

    sim = AcadosSim()

    # chain parameters
    chain_params = get_chain_params()
    chain_params["n_mass"] = args.n_masses
    n_mass = chain_params["n_mass"]
    M = chain_params["n_mass"] - 2 # number of intermediate masses
    # Ts = chain_params["Ts"]
    # Tsim = chain_params["Tsim"]
    N = chain_params["N"]
    # u_init = chain_params["u_init"]
    with_wall = chain_params["with_wall"]
    yPosWall = chain_params["yPosWall"]
    m = chain_params["m"]
    D = chain_params["D"]
    L = chain_params["L"]

    # export model
    model = export_chain_mass_model(n_mass, m, D, L)


    #%% simulate model
    sim.model = model
    sim.name = "chain_dirk_exp"

    Tf = 0.1
    nx = model.x.rows()
    nu = model.u.rows()
    N = 200

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.num_steps = 3
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.newton_iter = 5 # for implicit integrator>
    sim.solver_options.sens_forw = True
    sim.solver_options.sim_method_jac_reuse = args.jac_reuse
    sim.solver_options.exploit_semi_explicit = True
    sim.solver_options.num_stages = args.num_stages

    if COLLOC_TYPE == "GAUSS_RADAU_IIA":
        sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    elif COLLOC_TYPE == "SDIRK":
        sim.solver_options.collocation_type = "SDIRK"
        sim.solver_options.exploit_A_structure = args.exploit_dirk
    elif COLLOC_TYPE == "ESDIRK":
        sim.solver_options.collocation_type = "ESDIRK"
        sim.solver_options.exploit_A_structure = args.exploit_dirk

    # create
    acados_integrator = AcadosSimSolver(sim)

    simX = np.zeros((N+1, nx))

    # position of last mass
    xPosFirstMass = np.zeros((3,1))
    xEndRef = np.zeros((3,1))
    xEndRef[0] = L * (M+1) * 6
    # initial state
    pos0_x = np.linspace(xPosFirstMass[0], xEndRef[0], n_mass)
    x0 = np.zeros((nx, 1))
    x0[:3*(M+1):3] = pos0_x[1:].reshape((M+1,1))

    u0 = np.zeros((nu, 1))
    acados_integrator.set("u", u0)

    simX[0,:] = x0.flatten()

    t_total = []
    t_ad = []
    t_la = []
    t_factorize = []
    t_backsolve = []
    t_3 = []

    for i in range(N):
        # set initial state
        acados_integrator.set("x", simX[i,:])
        # solve
        status = acados_integrator.solve()
        # get solution
        simX[i+1,:] = acados_integrator.get("x")
        t_total.append(acados_integrator.get("time_tot"))
        t_ad.append(acados_integrator.get("time_ad"))
        t_la.append(acados_integrator.get("time_la"))
        t_factorize.append(acados_integrator.get("T1"))
        t_backsolve.append(acados_integrator.get("T2"))
        t_3.append(acados_integrator.get("T3"))


    print()
    print(f"colloc: {COLLOC_TYPE}, ns = {sim.solver_options.num_stages}, nx = {nx}")
    print(f"mean total time: {np.mean(t_total)}")
    print(f"mean ad time: {np.mean(t_ad)}")
    print(f"mean la time: {np.mean(t_la)}")
    print(f"mean factorize time: {np.mean(t_factorize)}")
    print(f"mean backsolve time: {np.mean(t_backsolve)}")
    print(f"mean t3 time: {np.mean(t_3)}")

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    if PLOT:
        plot_chain_position_traj(simX)

        #%% find steady state

        xrest = compute_steady_state(n_mass, m, D, L, xPosFirstMass, xEndRef)

        plot_chain_position(xrest, xPosFirstMass)

        plot_chain_position_3D((xrest, simX[-1,:]), xPosFirstMass)

        plt.show()

        #ani = animate_chain_position_3D(simX, xPosFirstMass)
        #ani = animate_chain_position(simX, xPosFirstMass)

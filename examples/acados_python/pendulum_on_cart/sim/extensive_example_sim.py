# -*- coding: future_fstrings -*-
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
sys.path.insert(0, '../common')

from acados_template import AcadosSim, AcadosSimSolver, acados_dae_model_json_dump, sim_get_default_cmake_builder
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np

sim = AcadosSim()

# model
model = export_pendulum_ode_model()
sim.model = model

Tf = 0.1
nx = model.x.size()[0]
nu = model.u.size()[0]
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

DETECT_GNSF = True
if sim.solver_options.integrator_type == "GNSF" and not DETECT_GNSF:
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

# create
cmake_builder = sim_get_default_cmake_builder()
acados_integrator = AcadosSimSolver(sim, cmake_builder=cmake_builder)

simX = np.ndarray((N+1, nx))
x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
u0 = np.array([0.0])
acados_integrator.set("u", u0)

simX[0,:] = x0


## Single test call
import time

t0 = time.time()
acados_integrator.set("seed_adj", np.ones((nx, 1)))
acados_integrator.set("x", x0)
acados_integrator.set("u", u0)
status = acados_integrator.solve()
time_external = time.time() - t0

S_forw = acados_integrator.get("S_forw")
Sx = acados_integrator.get("Sx")
Su = acados_integrator.get("Su")
S_hess = acados_integrator.get("S_hess")
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
print("S_hess, second order sensitivities:\n", S_hess)

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

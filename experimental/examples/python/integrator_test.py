#
#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren, Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor, Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan, Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
from acados import *
import time
import numpy as np
from models import pendulum_model, chen_model

# get model
ode_fun, nx, nu, impl_ode_fun = pendulum_model()

# define initial value, controls
x0 = np.array([ 4.2, 0.42, 0, 0])
u = np.array([1])

# create integrator
opts2 = {'step_size': 0.01}
start = time.time()
sim2 = integrator(ode_fun, opts2)
stop = time.time()
print("time to create" + str(stop-start))
# call integrator
print("integrator result: " + str(sim2.integrate(x0, u)))
# input("press any key to continue")

print("sim2 settings")
sim2.print_settings()

# input("press any key to continue")


# create another integrator
opts1 = {'step_size'        : 0.01, # only mandatory argument
        ### OPTIONAL ARGUMENTS
         'model_type'       : 0, # 0 - EXPLICIT (default)
                                 # 1 - IMPLICIT
         'integrator'       : "IRK", # default (ERK)
         'use_MX'           : False,  # default (False)
         'ns'               : 4, # default set in C implementation of integrator used
         'num_steps'        : 1, # default set in C implementation of integrator used
         'newton_iter'      : 3, # default set in C implementation of integrator used
         'output_z'         : 0, # default set in C implementation of integrator used
         'sens_forw'        : 1, # default set in C implementation of integrator used
         'sens_adj'         : 1, # default set in C implementation of integrator used
         'sens_hess'        : 1, # default set in C implementation of integrator used
         'sens_algebraic'   : 0, # default set in C implementation of integrator used
         'jac_reuse'        : 0, # default set in C implementation of integrator used
        }

sim1 = integrator(ode_fun, opts1)
print("sim1 created")

sim1.print_settings()

# set experiment parameters
M = 6
x1 = x0
x2 = x0
start_time = 0.0
t1 = 0.0
t2 = 0.0

# run experiment
for k in range(M):
    start_time = time.time()    # start timer
    x1 = np.array(sim1.integrate(x1, u))
    t1 +=  time.time() - start_time

    start_time = time.time()    # start timer
    x2 = np.array(sim2.integrate(x2, u))
    t2 +=  time.time() - start_time

#     sim1.set_step(sim1.step()*0.99)
    # step can be updated
    print("")
    print(x1)
    print(x2)

print("avg time IRK: ", t1/M)
print("avg time ERK: ", t2/M)

from acados import *
import time
import numpy as np
from models import pendulum_model, chen_model

# get model
ode_fun, nx, nu, impl_ode_fun = pendulum_model()

# create integrator
opts2 = {'step': 0.01}
sim2 = integrator(ode_fun, opts2)



# create another integrator
opts1 = {'step'             : 0.01, # only mandatory argument
        ### OPTIONAL ARGUMENTS
         'model_type'       : 1, # 0 - EXPLICIT (default)
                                 # 1 - IMPLICIT
         'integrator'       : "IRK", # default (ERK)
         'use_MX'           : True,  # default (False)
         'ns'               : 3, # default in integrator
         'num_steps'        : 1, # default in integrator
         'newton_iter'      : 3, # default in integrator
         'output_z'         : 0, # default in integrator
         'sens_forw'        : 0, # default in integrator
         'sens_adj'         : 0, # default in integrator
         'sens_hess'        : 0, # default in integrator
         'sens_algebraic'   : 0, # default in integrator
         'jac_reuse'        : 0, # default in integrator
        }
sim1 = integrator(impl_ode_fun, opts1)



# set initial value
x1 = np.array([ 4.2, 0.42, 0, 0])
x2 = x1
u = np.array([1])

# set experiment parameters
M = 100
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

    # sim1.set_step(sim1.step()*0.99)
    # step can be updated
    print("")
    print(x1)
    print(x2)

print("avg time sim1: ", t1/M)
print("avg time sim2: ", t2/M)

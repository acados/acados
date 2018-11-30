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

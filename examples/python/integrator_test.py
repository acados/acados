from acados import *
import time
import numpy as np
from models import pendulum_model, chen_model


ode_fun, nx, nu, impl_ode_fun = pendulum_model()

opts1 = {'step'      : 0.1, # only mandatory argument
         'sens_forw' : 0,   # OPTIONAL ARGUMENTS
         'sens_adj'  : 0,
         'sens_hess' : 0,
         'jac_reuse' : 0,
         'integrator': "IRK",
         'model_type': 1,
         'num_steps' : 1,
         'ns'        : 3}

sim1 = integrator(impl_ode_fun, opts1)

opts2 = {'step': 0.2, 'use_MX': True}
ode_fun, _, _ = chen_model()
sim2 = integrator(ode_fun, opts2)

x1 = np.array([ 4.2, 0.42, 0, 0])
u = np.array([1])

x2 = np.array([ 4.2, 0])

M = 2
start_time = 0.0
t1 = 0.0
t2 = 0.0
print("starting loop")
for k in range(M):
    start_time = time.time()    # start timer
    x1 = np.array(sim1.integrate(x1, u))
    t1 +=  time.time() - start_time
    print("CALLED integrator1 sucessfully")

    start_time = time.time()    # start timer
    x2 = np.array(sim2.integrate(x2, u))
    t2 +=  time.time() - start_time

    sim1.set_step(sim1.step()*0.99)

    print(x1)
    print(x2)

print("avg time sim1: ", t1/M)
print("avg time sim2: ", t2/M)

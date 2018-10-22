from acados import *
import time
import numpy as np
from models import pendulum_model, chen_model


ode_fun, nx, nu = pendulum_model()


opts1 = {'step'     : 0.1, # only mandatory argument
         'sens_forw': 1,
         'sens_adj' : 1,
         'jac_reuse': 0}
sim1 = integrator(ode_fun, opts1)

opts2 = {'step': 0.2, 'use_MX': True}
ode_fun, nx, nu = chen_model()
sim2 = integrator(ode_fun, opts2)

x1 = np.array([ 4.2, 0.42, 0, 0])
u = np.array([1])

x2 = np.array([ 4.2, 0])

M = 100
start_time = 0.0
t1 = 0.0
t2 = 0.0

for k in range(M):
    start_time = time.time()    # start timer
    x1 = np.array(sim1.integrate(x1, u))
    t1 +=  time.time() - start_time

    start_time = time.time()    # start timer
    x2 = np.array(sim2.integrate(x2, u))
    t2 +=  time.time() - start_time

    sim1.set_step(sim1.step()*0.99)

    print(x1)
    print(x2)

print("avg time sim1: ", t1/M)
print("avg time sim2: ", t2/M)

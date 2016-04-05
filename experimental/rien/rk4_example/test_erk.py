import erk_integrator as erk
import numpy as np

import matplotlib.pyplot as	plt

sim_in2 = erk.sim_in()
sim_out2 = erk.sim_out()

sim_in2.step = 0.1
sim_in2.nSteps = 1
Ts = 0.1
sim_in2.x = [1.0, -0.5]
sim_in2.u = [0.1]

print "x0 = ", sim_in2.x
print "u = ", sim_in2.u

# Perform a numerical simulation
N = 10;
xs2 = np.zeros((N+1,2))
xs2[0,:] = sim_in2.x
for i in range(N):
	sim_in2.x = xs2[i,:]
	erk.integrate(sim_in2, sim_out2)
	xs2[i+1,:] = sim_out2.xn

print "xnext = ", xs2[-1]

time = [Ts*i for i in range(N+1)]

plt.figure(1)
plt.clf()
plt.plot(time, xs2[:,0], 'o-', label="x1 (ERK)")
plt.plot(time, xs2[:,1], 'o-', label="x2 (ERK)")
plt.legend(loc='upper left')
plt.xlabel("time (s)")
plt.ylabel("state")
plt.show()

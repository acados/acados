import erk_integrator as erk
import sparse_erk_integrator as serk
import rk4_integrator as rk4
import numpy as np

import matplotlib.pyplot as	plt

sim_in = rk4.sim_in()
sim_out = rk4.sim_out()

sim_in2 = erk.sim_in()
sim_out2 = erk.sim_out()

sim_in3 = serk.sim_in()
sim_out3 = serk.sim_out()

sim_in.step = 0.1
sim_in2.step = 0.1
sim_in.nSteps = 1
sim_in2.nSteps = 1
Ts = sim_in.step*sim_in.nSteps

sim_in.x = [1.0, -0.5]
sim_in.u = [0.1]
sim_in2.x = [1.0, -0.5]
sim_in2.u = [0.1]

print "x0 = ", sim_in.x
print "u = ", sim_in.u

# Perform a numerical simulation
N = 10;
xs = np.zeros((N+1,2))
xs[0,:] = sim_in.x
xs2 = np.zeros((N+1,2))
xs2[0,:] = sim_in.x
Sxs = []
Sus = []
for i in range(N):
	sim_in.x = xs[i,:]
	rk4.integrate(sim_in, sim_out)
	xs[i+1,:] = sim_out.xn
	Sxs.append(sim_out.Sx)
	Sus.append(sim_out.Su)

	sim_in2.x = xs[i,:]
	erk.integrate(sim_in2, sim_out2)
	xs2[i+1,:] = sim_out2.xn

print "xnext = ", xs[-1]
print "Sx = ", Sxs[-1]
print "Su = ", Sus[-1]

time = [Ts*i for i in range(N+1)]

plt.figure(1)
plt.clf()
plt.plot(time, xs[:,0], 'o-', label="x1 (RK4)")
plt.plot(time, xs[:,1], 'o-', label="x2 (RK4)")
plt.plot(time, xs2[:,0], 'o-', label="x1 (ERK)")
plt.plot(time, xs2[:,1], 'o-', label="x2 (ERK)")
plt.legend(loc='upper	left')
plt.xlabel("time (s)")
plt.ylabel("state")
plt.show()


# Compute relative errors for varying step sizes
stepSize = 0.00001
numSteps = 100000
sim_in.step = stepSize
sim_in.nSteps = numSteps
sim_in.x = [1.0, -0.5]
sim_in.u = [0.1]
rk4.integrate(sim_in, sim_out);
x_ex = sim_out.xn  # 'exact' values

sim_in2.step = sim_in.step
sim_in2.nSteps = sim_in.nSteps
sim_in2.x = sim_in.x
sim_in2.u = sim_in.u

sim_in3.step = sim_in.step
sim_in3.nSteps = sim_in.nSteps
sim_in3.x = sim_in.x
sim_in3.u = sim_in.u

stepSize = 0.0001
numSteps = 10000
print "Ts = ", stepSize*numSteps

steps = [1, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
nSteps = [1, 2, 4, 10, 20, 40, 100, 200, 400, 1000]

xs3 = np.zeros((len(steps),2))
err = np.zeros(len(steps))
timings = np.zeros(len(steps))
timings2 = np.zeros(len(steps))
timings3 = np.zeros(len(steps))
for i in range(len(steps)):
    sim_in.step = steps[i]
    sim_in.nSteps = nSteps[i]
    rk4.integrate(sim_in, sim_out);
    xs3[i,:] = sim_out.xn
    err[i] = max(abs(xs3[i,:]-x_ex))
    timings[i] = 1e6*sim_out.cpuTime

    sim_in2.step = steps[i]
    sim_in2.nSteps = nSteps[i]
    erk.integrate(sim_in2, sim_out2);
    timings2[i] = 1e6*sim_out2.cpuTime

    sim_in3.step = steps[i]
    sim_in3.nSteps = nSteps[i]
    serk.integrate(sim_in3, sim_out3);
    timings3[i] = 1e6*sim_out3.cpuTime

print "abs err = ", err

plt.figure(2)
plt.clf()
plt.plot(nSteps, err[:], 'ro--')
plt.xlabel("number of steps")
plt.ylabel("max abs error")
plt.yscale('log')
plt.xscale('log')
plt.show()


plt.figure(3)
plt.clf()
plt.plot(nSteps[1:8], timings[1:8], 'go-', label="RK4 implementation")
plt.plot(nSteps[1:8], timings2[1:8], 'ro-', label="ERK implementation")
plt.plot(nSteps[1:8], timings3[1:8], 'bo-', label="'SPARSE' ERK implementation")
plt.legend(loc='upper	left')
plt.xlabel("number of steps")
plt.ylabel("cpuTime/step ($\mu$s)")
plt.show()

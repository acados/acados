import matplotlib.pyplot as plt
from numpy import array, concatenate, cos, linspace, pi, reshape

from acados import integrator
from models import pendulum_model

ode_fun, nx, _ = pendulum_model()

time_grid, dt = linspace(start=0, stop=10, num=51, retstep=True)
sim_options = {'time_step': dt, 'order': 4}

integrator = integrator('erk', ode_fun, sim_options)

current_state = array([0.0, pi, 0.0, 0.0]) # pendulum hangs down
simulation = [current_state]

for t in time_grid[:-1]:
    # apply periodic force to pendulum with period 5s.
    output = integrator.evaluate(current_state, cos(2*pi*t/5))
    current_state = output.final_state
    simulation.append(current_state)

simulation = reshape(concatenate(simulation), (-1, nx))

plt.ion()
plt.plot(time_grid, simulation)
plt.legend(['p', 'theta', 'v', 'omega'])
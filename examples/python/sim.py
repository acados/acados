from numpy import array

from acados import sim_solver
from example_model import example_model

_, _, ode_fun = example_model()

sim_options = {'time_step': 0.2, 'order': 4}

sim = sim_solver('rk', ode_fun, sim_options)

current_state = [1.0, 1.0]
control = 0.25
output = sim.evaluate(current_state, control)

assert(all(abs(output.final_state - array([1.2670, 1.1442])) < 1e-3))
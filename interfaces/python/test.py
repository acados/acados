from casadi import *
import acados_integrator as ai
import time


nx = 4
nu = 1

x = SX.sym('x', nx, 1)
ode_expr = -2*x


start_time = time.time()    # start timer

sim_model = ai.acados_integrator_model()
sim_model.set('type', 'explicit')
sim_model.set('ode_expr', ode_expr)
sim_model.set('x', x)
sim_model.set('model_name', 'model')

end_time = time.time()      # stop timer
print('sim_model time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

sim_opts = ai.acados_integrator_opts()
sim_opts.set('scheme', 'erk')
sim_opts.set('sens_forw', 'false')

end_time = time.time()      # stop timer
print('sim_opts time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

sim = ai.acados_integrator(sim_model, sim_opts)

end_time = time.time()      # stop timer
print('sim create time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

x0 = np.array([1, 0, 2, -1])
sim.set('x', x0)
sim.set('t', 0.05)

end_time = time.time()      # stop timer
print('sim set x time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

flag = sim.solve()

end_time = time.time()      # stop timer
print('sim solve time {:e}'.format(end_time-start_time))
print(flag)



start_time = time.time()    # start timer

xn = sim.get('xn')

end_time = time.time()      # stop timer
print('sim get xn time {:e}'.format(end_time-start_time))
print(xn)

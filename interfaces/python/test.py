from casadi import *
import acados_integrator as ai
import time


nx = 4
#nu = 0

x = SX.sym('x', nx, 1)
xdot = SX.sym('xdot', nx, 1)

expl_ode_expr = -2*x
impl_ode_expr = xdot - expl_ode_expr


start_time = time.time()    # start timer

sim_model = ai.acados_integrator_model()
#sim_model.set('type', 'explicit')
sim_model.set('type', 'implicit')
#sim_model.set('ode_expr', expl_ode_expr)
sim_model.set('ode_expr', impl_ode_expr)
sim_model.set('x', x)
sim_model.set('xdot', xdot)
#sim_model.set('model_name', 'expl_model')
sim_model.set('model_name', 'impl_model')

end_time = time.time()      # stop timer
print('sim_model time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

sim_opts = ai.acados_integrator_opts()
#sim_opts.set('scheme', 'erk')
sim_opts.set('scheme', 'irk')
sim_opts.set('sens_forw', 'true')
#sim_opts.set('sens_forw', 'false')

end_time = time.time()      # stop timer
print('sim_opts time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

sim = ai.acados_integrator(sim_model, sim_opts)

end_time = time.time()      # stop timer
print('sim create time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

x0 = np.array([1, 0, 2, -1])
xdot0 = x0 + 2*x0
sim.set('x', x0)
sim.set('xdot', xdot0)
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
Sxn = sim.get('Sxn')

end_time = time.time()      # stop timer
print('sim get xn time {:e}'.format(end_time-start_time))
print(xn)
print(Sxn)

from casadi import *
from acados.sim import *
#import pickle as pkl
import time





## create one integrator

print('\nfirst integrator\n')

nx = 4
#nu = 0

x = SX.sym('x', nx, 1)
xdot = SX.sym('xdot', nx, 1)

expl_ode_expr = -2*x
impl_ode_expr = xdot - expl_ode_expr

start_time = time.time()    # start timer

sim_model = acados_integrator_model()
#sim_model.set('model_name', 'expl_model')
sim_model.set('model_name', 'impl_model')
#sim_model.set('type', 'explicit')
sim_model.set('type', 'implicit')
#sim_model.set('ode_expr', expl_ode_expr)
sim_model.set('ode_expr', impl_ode_expr)
sim_model.set('x', x)
sim_model.set('xdot', xdot)

end_time = time.time()      # stop timer
print('sim_model time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

sim_opts = acados_integrator_opts()
#sim_opts.set('scheme', 'erk')
sim_opts.set('scheme', 'irk')
#sim_opts.set('sens_forw', 'true')
sim_opts.set('sens_forw', 'false')
sim_opts.set('codgen_model', 'true')
#sim_opts.set('codgen_model', 'false')

end_time = time.time()      # stop timer
print('sim_opts time {:e}'.format(end_time-start_time))



start_time = time.time()    # start timer

sim = acados_integrator(sim_model, sim_opts)

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





## create a sendond integrator

print('\nsecond integrator\n')

# update expression
expl_ode_expr = -1*x
impl_ode_expr = xdot - expl_ode_expr

# update model
#sim_model.set('model_name', 'new_impl_model')
#sim_model.set('type', 'implicit')
sim_model.set('ode_expr', impl_ode_expr)

# update opts
#sim_opts.set('scheme', 'irk')
#sim_opts.set('sens_forw', 'true')
#sim_opts.set('codgen_model', 'true')

# new sim solver
new_sim = acados_integrator(sim_model, sim_opts)

new_sim.set('x', x0)
new_sim.set('xdot', xdot0)
new_sim.set('t', 0.05)

flag = new_sim.solve()

xn = new_sim.get('xn')
Sxn = new_sim.get('Sxn')

print(xn)
print(Sxn)









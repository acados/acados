from casadi import *
import acados_integrator as ai


nx = 4
nu = 1

x = SX.sym('x', nx, 1)
ode_expr = -2*x


sim_model = ai.acados_integrator_model()
sim_model.set('ode_expr', ode_expr)
sim_model.set('x', x)
sim_model.set('model_name', 'model')

#sim_model.generate_lib('model')



sim = ai.acados_integrator(sim_model)

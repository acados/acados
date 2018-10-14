from casadi import *
import acados_integrator as ai


nx = 4
nu = 1

x = SX.sym('x', nx, 1)
casadi_ode_expr = -2*x

# Form a function
user_fun_name = 'ode_expr'
python_ode_expr = Function(user_fun_name, [x], [casadi_ode_expr], ['x'], ['ode_expr_out'])



sim_model = ai.acados_integrator_model()
sim_model.set('ode_expr', python_ode_expr)
sim_model.set('x', x)

sim_model.generate_lib('model')



sim = ai.acados_integrator(sim_model)

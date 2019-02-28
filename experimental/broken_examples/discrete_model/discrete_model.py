from casadi import *

x = SX.sym('x')
u = SX.sym('u')
xu = vertcat(x, u)

rhs = x - u

discrete_model = Function('discrete_model', [x, u], [rhs, jacobian(rhs, xu)])
discrete_model.generate('discrete_model.c', {'with_header': True})

discrete_model_cost = Function('discrete_model_cost', [x, u], [xu, jacobian(xu, xu)])
discrete_model_cost.generate('discrete_model_cost.c', {'with_header': True})

discrete_model_costN = Function('discrete_model_costN', [x], [x, jacobian(x, x)])
discrete_model_costN.generate('discrete_model_costN.c', {'with_header': True})
import matplotlib.pyplot as plt
from numpy import array, diag
from scipy.linalg import block_diag

from acados import *
from casadi import *

x = SX.sym('x',2);
u = SX.sym('u',2);
y = vertcat(mtimes(x.T,x) + mtimes(u.T,u), 3*mtimes(u.T,u)) ;
yfun = Function('y',[x,u],[y]);

nlp_function = ocp_nlp_function(yfun);

x_in = [1.0,2.0];
u_in = [2.0,3.0];

res = nlp_function.evaluate(x_in, u_in);

print(res)

print(ls_cost.ls_cost_matrix)
print(ls_cost.ls_cost_ref)

from acados import *
from numpy import array, ones, zeros, diag

qp_in = ocp_qp({'N':5, 'nx':2, 'nu':1})
x0 = array([1,1])

# specify OCP
qp_in.A = array([[0,1],[0,0]])
qp_in.B = array([[0],[1]])
qp_in.Q = diag([1,1])
qp_in.R = diag([1])

# specify initial condition
x0 = array([1,1])
qp_in.lb[0] = x0
qp_in.ub[0] = x0

solver = ocp_qp_solver("condensing_qpoases", qp_in)
result = solver.solve()
print(result)

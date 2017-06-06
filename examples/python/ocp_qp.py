from acados import *
from numpy import array, ones, zeros, diag

qp = ocp_qp({'N':5, 'nx':2, 'nu':1})

# specify OCP
qp.A = array([[0,1],[0,0]])
qp.B = array([[0],[1]])
qp.Q = diag([1,1])
qp.R = diag([1])

# specify initial condition
x0 = array([1,1])
qp.lb[0] = x0
qp.ub[0] = x0

# solve QP
solver = ocp_qp_solver("condensing_qpoases", qp)
result = solver.solve()
assert(abs(-0.5 - result[0]) < 1e-8)
print(result)

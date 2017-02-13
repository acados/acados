from acados import *
from numpy import array, ones, zeros, diag

N = 5
nx = 2
nu = 1
nb = [nx] + N*[0]

x0 = array([1,1])

qp_in = ocp_qp_in({'N':5, 'nx':nx, 'nu':nu, 'nb':nb})
qp_in.idxb[0] = array([1, 2])

qp_in.lb[0] = x0
qp_in.ub[0] = x0

A = zeros([nx, nx])
A[0,1] = 1
qp_in.A = A
B = zeros([nx, nu])
B[1,0] = 1
qp_in.B = B

qp_in.Q = diag([1,1])
qp_in.R = diag([1])

solver = ocp_qp_solver("condensing_qpoases", qp_in)

result = solver.solve()

print(result)

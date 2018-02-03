from numpy import array, diag

from acados import ocp_qp

qp = ocp_qp(N=5, nx=2, nu=1)

# specify OCP
qp.update('A', array([[0, 1], [0, 0]]))
qp.update('B', array([[0], [1]]))
qp.update('Q', diag([1, 1]))
qp.update('R', diag([1]))

# specify initial condition
x0 = array([1, 1])
qp.update('lbx', 0, x0)
qp.update('ubx', 0, x0)

# solve QP
output = qp.solve("qpoases")
# Inspect the output struct
print(output.states())
print(output.info())

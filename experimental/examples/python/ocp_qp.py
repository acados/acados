
from numpy import array, diag, inf

from acados import ocp_qp

qp = ocp_qp(N=5, nx=2, nu=1)

# specify OCP
qp.set_field('A', array([[1, 1], [0, 1]]))
qp.set_field('B', array([[0], [1]]))
qp.set_field('Q', diag([1, 1]))
qp.set_field('R', diag([1]))

# specify bounds
qp.set_field('lbx', array([0.5, -inf]))
qp.set_field('ubx', array([3.0, +inf]))

# specify initial condition
x0 = array([1.1, 1.1])
qp.set_field('lbx', 0, x0)
qp.set_field('ubx', 0, x0)

# solve QP
qp.initialize_solver('qpoases')
output = qp.solve()

print("Optimal state trajectory:")
print(output.states())

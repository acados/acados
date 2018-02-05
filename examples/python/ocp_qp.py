from numpy import array, diag, inf

from acados import ocp_qp

qp = ocp_qp(N=5, nx=2, nu=1, nbx=2)

# specify OCP
qp.set('A', array([[0, 1], [0, 0]]))
qp.set('B', array([[0], [1]]))
qp.set('Q', diag([1, 1]))
qp.set('R', diag([1]))

# specify bounds
qp.set("lbx", array([0.5, 0.5]))
qp.set("ubx", array([2, 2]))

# specify initial condition
x0 = array([1, 1])
qp.set('lbx', 0, x0)
qp.set('ubx', 0, x0)

# solve QP
output = qp.solve("qpoases", {'max_nwsr': 5})
# Inspect the output struct
print(output.states())
print(output.info())

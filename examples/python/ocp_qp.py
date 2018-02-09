
from numpy import array, diag, inf

from acados import ocp_qp

qp = ocp_qp(N=5, nx=2, nu=1)

# specify OCP
qp.set('A', array([[0, 1], [0, 0]]))
qp.set('B', array([[0], [1]]))
qp.set('Q', diag([1, 1]))
qp.set('R', diag([1]))

# specify bounds
# qp.set("lbx", array([0.5, -inf]))
# qp.set("ubx", array([2.0, +inf]))

# qp.set("lbu", array([-inf]))
# qp.set("ubu", array([+inf]))

# specify initial condition
x0 = array([1, 1])
qp.set('lbx', 0, x0)
qp.set('ubx', 0, x0)

# solve QP and print solution
# for solver_name in ("sparse_hpipm", "condensing_hpipm", "hpmpc", "qpdunes", "qore", "qpoases"):
for solver_name in ("sparse_hpipm",):
    print(solver_name + ": ")
    qp.initialize_solver(solver_name)
    output = qp.solve()
    print(output.states())
    print(output.info())
    print()

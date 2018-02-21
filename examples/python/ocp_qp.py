
from numpy import array, diag, inf

from acados import ocp_qp

qp = ocp_qp(N=5, nx=2, nu=1)

# specify OCP
qp.set('A', array([[1, 1], [0, 1]]))
qp.set('B', array([[0], [1]]))
qp.set('Q', diag([1, 1]))
qp.set('R', diag([1]))

# specify bounds
# qp.set("lbx", array([0.5, -inf]))
# qp.set("ubx", array([2.0, +inf]))

qp.set("q", array([1.0, 1.0]))
# specify initial condition
x0 = array([1, 1])
qp.set('lbx', 0, x0)
qp.set('ubx', 0, x0)

qp.initialize_solver("sparse_hpipm")
# solve QP and print solution
output = qp.solve()


qp.set("lbx", 1, array([0.5, 0.5]))
qp.set("ubx", 1, array([2.0, 2.0]))

qp.initialize_solver("sparse_hpipm")

output = qp.solve()


# for solver_name in ("sparse_hpipm", "condensing_hpipm", "hpmpc", "qpdunes", "qore", "qpoases"):
    # print(solver_name + ": ")
    # qp.set("lbx", 1, array([0.5, 2.0]))
    # qp.set("ubx", 1, array([0.5, 2.0]))
    # output = qp.solve()
    # print(output.states())
    # print(output.info())
    # print()


from numpy import array, diag, inf

from acados import ocp_qp

def solve(qp):
    for solver_name in ("qpdunes",):
        print(solver_name + ": ")
        qp.initialize_solver(solver_name)
        output = qp.solve()
        print(output.states())
        print(output.controls())
        print(output.info())
        print()

qp = ocp_qp(N=5, nx=2, nu=1)

# specify OCP
qp.set_field('A', array([[1, 1], [0, 1]]))
qp.set_field('B', array([[0], [1]]))
qp.set_field('Q', diag([1, 1]))
qp.set_field('R', diag([1]))

# specify bounds
# qp.set_field("lbx", array([0.5, -inf]))
# qp.set_field("ubx", array([2.0, +inf]))

# qp.set_field("q", array([1.0, 1.0]))
# specify initial condition
x0 = array([1.1, 1.1])
qp.set_field('lbx', 0, x0)
qp.set_field('ubx', 0, x0)

# solve(qp)

# qp.set_field("lbx", 3, array([1.0, 1.0]))
# qp.set_field("ubx", 3, array([2.0, +inf]))

# solve(qp)

# qp.set_field("lbx", 2, array([-inf, -10.0]))
# qp.set_field("ubx", 2, array([+inf, +10.0]))
qp.initialize_solver("sparse_hpipm")
output = qp.solve()
print(output.states())

solve(qp)

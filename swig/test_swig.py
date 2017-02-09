import numpy

print('Numpy version:')
print(numpy.version.full_version)

from acados import *

nx = 4
nu = 2
nb = 1

qp_in = ocp_qp_in({'N':5, 'nx':nx, 'nu':nu, 'nb':nb})
qp_out = ocp_qp_out()
qp_arg = ocp_qp_condensing_qpoases_args()
qp_mem = None

# initialise_qpoases(qp_in)
# ocp_qp_condensing_qpoases(qp_in, qp_out, qp_arg, qp_mem)

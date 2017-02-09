from acados import *
import numpy

nx = 4
nu = 2
nb = 1

qp = ocp_qp_in({'N':5, 'nx':nx, 'nu':nu, 'nb':nb})

print(qp.r)

qp.r = numpy.ones(nu)

print(qp.r)

print(qp.A)

qp.A = numpy.ones([nx,nx])

print(qp.r)

qp.r[0] = numpy.zeros(nu)

print(qp.r)

print(qp.idxb)

qp.idxb = numpy.ones(nb)

print(qp.idxb)

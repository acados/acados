from ctypes import *

acados   = CDLL('c_generated_code/acados_solver_pendulum_ode.so')
import pdb; pdb.set_trace()
acados.acados_create()
acados.call_solve()
acados.call_destroy()

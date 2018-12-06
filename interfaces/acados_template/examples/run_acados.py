from ctypes import *

acados   = CDLL('c_generated_code/libacados_pendulum.so')
acados.call_acados()

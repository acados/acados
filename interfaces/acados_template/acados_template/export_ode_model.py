from casadi import *
class ode_model():
    def __init__(self):
        self.f_impl_expr = None
        self.f_expl_expr = None
        self.x = None
        self.xdot = None
        self.u = None
        self.z = None
        self.name = None


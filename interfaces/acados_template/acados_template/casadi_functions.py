from casadi import *
class acados_dae():
    def __init__(self):
        self.f_impl_expr = None
        self.f_expl_expr = None
        self.x = None
        self.xdot = None
        self.u = None
        self.z = None
        self.name = None

class acados_constraint():
    def __init__(self):
        self.expr = None
        self.x = None
        self.u = None
        self.z = None
        self.name = None

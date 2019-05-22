from casadi import *
class acados_dae():
    def __init__(self):
        self.f_impl_expr = None #: CasADi expression for the implicit dynamics :math:`F(\dot{x}, x, u, z) = 0`
        self.f_expl_expr = None #: CasADi expression for the explicit dynamics :math:`\dot{x} = f(x, u)`
        self.x = None           #: CasADi variable describing the state of the system
        self.xdot = None        #: CasADi variable describing the derivative of the state wrt time
        self.u = None           #: CasADi variable describing the input of the system
        self.z = None           #: CasADi variable describing the algebraic variables of the DAE
        self.name = None        #: name associated with the function

class acados_constraint():
    def __init__(self):
        self.expr = None #: CasADi expression for the constraint
        self.x = None    #: CasADi variable describing the state of the system
        self.u = None    #: CasADi variable describing the input of the system
        self.z = None    #: CasADi variable describing the algebraic variables of the DAE
        self.name = None #: name associated with the function


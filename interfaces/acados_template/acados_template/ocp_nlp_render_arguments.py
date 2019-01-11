import numpy as np

class ocp_nlp_dims:
    def __init__(self):
        self._nx  = None  # number of states
        self._nz  = 0     # number of algebraic variables
        self._nu  = None  # number of inputs
        self._ny  = None  # number of residuals in Lagrange term
        self._nyN = None  # number of residuals in Mayer term
        self._nbx = 0     # number of state bounds 
        self._nbu = 0     # number of input bounds
        self._nu  = None  # number of inputs
        self._N   = None  # prediction horizon 

    @property
    def nx(self):
        return self._nx

    @property
    def nz(self):
        return self._nz

    @property
    def nu(self):
        return self._nu

    @property
    def ny(self):
        return self._ny

    @property
    def nyN(self):
        return self._nyN

    @property
    def nbx(self):
        return self._nbx

    @property
    def nbu(self):
        return self._nbu
    
    @property
    def N(self):
        return self._N

    @nx.setter
    def nx(self, nx):
        if type(nx) == int and nx > 0:
            self._nx = nx
        else:
            raise Exception('Invalid nx value. Exiting.')

    @nz.setter
    def nz(self, nz):
        if type(nz) == int and nz > -1:
            self._nz = nz
        else:
            raise Exception('Invalid nz value. Exiting.')

    @ny.setter
    def ny(self, ny):
        if type(ny) == int and ny > 0:
            self._ny = ny
        else:
            raise Exception('Invalid ny value. Exiting.')

    @nyN.setter
    def nyN(self, nyN):
        if type(nyN) == int and nyN > 0:
            self._nyN = nyN
        else:
            raise Exception('Invalid nyN value. Exiting.')

    @nu.setter
    def nu(self, nu):
        if type(nu) == int and nu > 0:
            self._nu = nu
        else:
            raise Exception('Invalid nu value. Exiting.')

    @nbu.setter
    def nbu(self, nbu):
        if type(nbu) == int and nbu > -1:
            self._nbu = nbu
        else:
            raise Exception('Invalid nbu value. Exiting.')

    @nbx.setter
    def nbx(self, nbx):
        if type(nbx) == int and nbx > -1:
            self._nbx = nbx
        else:
            raise Exception('Invalid nbx value. Exiting.')

    @N.setter
    def N(self, N):
        if type(N) == int and N > 0:
            self._N = N
        else:
            raise Exception('Invalid N value. Exiting.')

class ocp_nlp_cost:
    # linear least-squares cost: || Vx*x + Vu*x + Vz*z ||^2_W
    def __init__(self):
        # Lagrange term
        self._W     = None  # weight matrix
        self._Vx    = None  # x matrix coefficient
        self._Vu    = None  # u matrix coefficient
        self._Vz    = None  # z matrix coefficient
        self._yref  = None  # reference
        # Mayer term
        self._WN    = None  # weight matrix
        self._VxN   = None  # x matrix coefficient
        self._yrefN = None  # reference

    # Lagrange term
    @property
    def W(self):
        return self._W

    @property
    def Vx(self):
        return self._Vx

    @property
    def Vu(self):
        return self._Vu

    @property
    def Vz(self):
        return self._Vz

    @property
    def yref(self):
        return self._yref

    @W.setter
    def W(self, W):
        if type(W) == np.ndarray:
            self._W = W
        else:
            raise Exception('Invalid W value. Exiting.')
    
    @Vx.setter
    def Vx(self, Vx):
        if type(Vx) == np.ndarray:
            self._Vx = Vx
        else:
            raise Exception('Invalid Vx value. Exiting.')
    
    @Vu.setter
    def Vu(self, Vu):
        if type(Vu) == np.ndarray:
            self._Vu = Vu
        else:
            raise Exception('Invalid Vu value. Exiting.')

    @Vz.setter
    def Vz(self, Vz):
        if type(Vz) == np.ndarray:
            self._Vz = Vz
        else:
            raise Exception('Invalid W value. Exiting.')

    @yref.setter
    def yref(self, yref):
        if type(yref) == np.ndarray:
            self._yref = yref
        else:
            raise Exception('Invalid yref value. Exiting.')

    # Mayer term
    @property
    def WN(self):
        return self._WN

    @property
    def VxN(self):
        return self._VxN

    @property
    def yrefN(self):
        return self._yrefN

    @WN.setter
    def WN(self, WN):
        if type(WN) == np.ndarray:
            self._WN = WN
        else:
            raise Exception('Invalid WN value. Exiting.')
    
    @VxN.setter
    def VxN(self, VxN):
        if type(VxN) == np.ndarray:
            self._VxN = VxN
        else:
            raise Exception('Invalid VxN value. Exiting.')

    @yrefN.setter
    def yrefN(self, yrefN):
        if type(yrefN) == np.ndarray:
            self._yrefN = yrefN
        else:
            raise Exception('Invalid yrefN value. Exiting.')
    
class ocp_nlp_constraints:
    def __init__(self):
        self._lbx = None  
        self._lbu = None  
        self._ubx = None  
        self._ubu = None  
        self._x0 = None  

    @property
    def lbx(self):
        return self._lbx

    @property
    def lbu(self):
        return self._lbu
    
    @property
    def ubx(self):
        return self._ubx

    @property
    def ubu(self):
        return self._ubu

    @property
    def x0(self):
        return self._x0

    @lbx.setter
    def lbx(self, lbx):
        if type(lbx) == np.ndarray:
            self._lbx = lbx
        else:
            raise Exception('Invalid lbx value. Exiting.')

    @lbu.setter
    def lbu(self, lbu):
        if type(lbu) == np.ndarray:
            self._lbu = lbu
        else:
            raise Exception('Invalid lbu value. Exiting.')

    @ubx.setter
    def ubx(self, ubx):
        if type(ubx) == np.ndarray:
            self._ubx = ubx
        else:
            raise Exception('Invalid ubx value. Exiting.')

    @ubu.setter
    def ubu(self, ubu):
        if type(ubu) == np.ndarray:
            self._ubu = ubu
        else:
            raise Exception('Invalid ubu value. Exiting.')

    @x0.setter
    def x0(self, x0):
        if type(x0) == np.ndarray:
            self._x0 = x0
        else:
            raise Exception('Invalid x0 value. Exiting.')

class ocp_nlp_solver_config:
    def __init__(self):
        self._qp_solver      = 'PARTIAL_CONDENSING_HPIPM'   # qp solver to be used in the NLP solver
        self._hessian_approx = 'GAUSS_NEWTON'               # hessian approximation
        self._integrator_type = 'ERK'                       # integrator type
        self._tf = None                                     # prediction horizon
        self._nlp_solver_tpye = 'SQP_RTI'                   # NLP solver 

    @property
    def qp_solver(self):
        return self._qp_solver

    @property
    def hessian_approx(self):
        return self._hessian_approx

    @property
    def integrator_type(self):
        return self._integrator_type

    @property
    def nlp_solver_type(self):
        return self._nlp_solver_type

    @qp_solver.setter
    def qp_solver(self, qp_solver):
        qp_solvers = ('PARTIAL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_QPOASES', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_QPOASES')

        if type(qp_solver) == str and qp_solver in qp_solvers:
            self._qp_solver = qp_solver
        else:
            raise Exception('Invalid qp_solver value. Possible values are:\n\n' \
                    + ',\n'.join(qp_solvers) + '.\n\nYou have: ' + qp_solver + '.\n\nExiting.')
    @property
    def tf(self):
        return self._tf

    @hessian_approx.setter
    def hessian_approx(self, hessian_approx):
        hessian_approxs = ('GAUSS_NEWTON', 'EXACT')

        if type(hessian_approx) == str and hessian_approx in hessian_approxs:
            self._hessian_approx = hessian_approx
        else:
            raise Exception('Invalid hessian_approx value. Possible values are:\n\n' \
                    + ',\n'.join(hessian_approxs) + '.\n\nYou have: ' + hessian_approx + '.\n\nExiting.')

    @integrator_type.setter
    def integrator_type(self, integrator_type):
        integrator_types = ('ERK', 'IRK')

        if type(integrator_type) == str and integrator_type in integrator_types:
            self._integrator_type = integrator_type
        else:
            raise Exception('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(integrator_types) + '.\n\nYou have: ' + integrator_type + '.\n\nExiting.')

    @tf.setter
    def tf(self, tf):
        self._tf = tf

    @nlp_solver_type.setter
    def nlp_solver_type(self, nlp_solver_type):
        nlp_solver_types = ('SQP', 'SQP_RTI')

        if type(nlp_solver_type) == str and nlp_solver_type in nlp_solver_types:
            self._nlp_solver_type = nlp_solver_type
        else:
            raise Exception('Invalid nlp_solver_type value. Possible values are:\n\n' \
                    + ',\n'.join(nlp_solver_types) + '.\n\nYou have: ' + nlp_solver_type + '.\n\nExiting.')

class ocp_nlp_constant:
    def __init__(self):
        self.name  = None # constant name
        self.value = None # constant value

class ocp_nlp_render_arguments:
    def __init__(self):
        self.dims = ocp_nlp_dims()
        self.cost = ocp_nlp_cost()
        self.constraints = ocp_nlp_constraints()
        self.solver_config = ocp_nlp_solver_config()
        self.model_name = None 
        self.constants = []
        self.acados_include_path = []
        self.acados_lib_path = []

def check_ra(ra):
    if ra.solver_config.hessian_approx == 'EXACT' and ra.solver_config.integrator_type == 'IRK':
        raise Exception('Exact Hessians not yet supported with IRK integrators.')

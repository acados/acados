import numpy as np

class ocp_nlp_dims:
    def __init__(self):
        self._nx   = None  # number of states
        self._nz   = 0     # number of algebraic variables
        self._nu   = None  # number of inputs
        self._np   = 0     # number of parameters
        self._ny   = None  # number of residuals in Lagrange term
        self._nyN  = None  # number of residuals in Mayer term
        self._npd  = 0     # number of positive definite constraints
        self._npdN = 0     # number of positive definite constraints in last stage
        self._nh   = 0     # number of nonlinear constraints
        self._nhN  = 0     # number of nonlinear constraints in last stage
        self._nbx  = 0     # number of state bounds 
        self._nbu  = 0     # number of input bounds
        self._ng   = 0     # number of general constraints
        self._nbxN = 0     # number of state bounds in last stage 
        self._ngN  = 0     # number of general constraints in last stage
        self._N    = None  # prediction horizon 

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
    def np(self):
        return self._np

    @property
    def ny(self):
        return self._ny

    @property
    def npd(self):
        return self._npd

    @property
    def npdN(self):
        return self._npdN

    @property
    def nh(self):
        return self._nh

    @property
    def nhN(self):
        return self._nhN

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
    def ng(self):
        return self._ng

    @property
    def nbxN(self):
        return self._nbxN
    
    @property
    def ngN(self):
        return self._ngN

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

    @nu.setter
    def nu(self, nu):
        if type(nu) == int and nu > 0:
            self._nu = nu
        else:
            raise Exception('Invalid nu value. Exiting.')

    @np.setter
    def np(self, np):
        if type(np) == int and np > -1:
            self._np = np
        else:
            raise Exception('Invalid np value. Exiting.')

    @ny.setter
    def ny(self, ny):
        if type(ny) == int and ny > -1:
            self._ny = ny
        else:
            raise Exception('Invalid ny value. Exiting.')

    @nyN.setter
    def nyN(self, nyN):
        if type(nyN) == int and nyN > -1:
            self._nyN = nyN
        else:
            raise Exception('Invalid nyN value. Exiting.')

    @npd.setter
    def npd(self, npd):
        if type(npd) == int and npd > -1:
            self._npd = npd
        else:
            raise Exception('Invalid npd value. Exiting.')

    @npdN.setter
    def npdN(self, npdN):
        if type(npdN) == int and npdN > -1:
            self._npdN = npdN
        else:
            raise Exception('Invalid npdN value. Exiting.')

    @nh.setter
    def nh(self, nh):
        if type(nh) == int and nh > -1:
            self._nh = nh
        else:
            raise Exception('Invalid nh value. Exiting.')

    @nhN.setter
    def nhN(self, nhN):
        if type(nhN) == int and nhN > -1:
            self._nhN = nhN
        else:
            raise Exception('Invalid nhN value. Exiting.')

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

    @ng.setter
    def ng(self, ng):
        if type(ng) == int and ng > -1:
            self._ng = ng
        else:
            raise Exception('Invalid ng value. Exiting.')

    @nbxN.setter
    def nbxN(self, nbxN):
        if type(nbxN) == int and nbxN > -1:
            self._nbxN = nbxN
        else:
            raise Exception('Invalid nbxN value. Exiting.')

    @ngN.setter
    def ngN(self, ngN):
        if type(ngN) == int and ngN > -1:
            self._ngN = ngN
        else:
            raise Exception('Invalid ngN value. Exiting.')

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
        self._lbx    = None  
        self._lbu    = None  
        self._idxbx  = None
        self._ubx    = None  
        self._ubu    = None  
        self._idxbu  = None
        self._lg     = None  
        self._ug     = None  
        self._lh     = None  
        self._uh     = None  
        self._D      = None  
        self._C      = None  
        self._lbxN   = None  
        self._ubxN   = None  
        self._idxbxN = None
        self._CN     = None  
        self._lgN    = None  
        self._ugN    = None  
        self._lhN    = None  
        self._uhN    = None  
        self._x0     = None  

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
    def idxbx(self):
        return self._idxbx

    @property
    def idxbu(self):
        return self._idxbu

    @property
    def lg(self):
        return self._lg

    @property
    def ug(self):
        return self._ug

    @property
    def lh(self):
        return self._lh

    @property
    def uh(self):
        return self._uh

    @property
    def D(self):
        return self._D

    @property
    def C(self):
        return self._C

    @property
    def lbxN(self):
        return self._lbxN

    @property
    def ubxN(self):
        return self._ubxN

    @property
    def idxbxN(self):
        return self._idxbxN

    @property
    def CN(self):
        return self._CN

    @property
    def lgN(self):
        return self._lgN

    @property
    def ugN(self):
        return self._ugN

    @property
    def lgN(self):
        return self._lgN

    @property
    def ugN(self):
        return self._ugN

    @property
    def x0(self):
        return self._x0

    @property
    def p(self):
        return self._p

    @lbx.setter
    def lbx(self, lbx):
        if type(lbx) == np.ndarray:
            self._lbx = lbx
        else:
            raise Exception('Invalid lbx value. Exiting.')

    @ubx.setter
    def ubx(self, ubx):
        if type(ubx) == np.ndarray:
            self._ubx = ubx
        else:
            raise Exception('Invalid ubx value. Exiting.')

    @idxbx.setter
    def idxbx(self, idxbx):
        if type(idxbx) == np.ndarray:
            self._idxbx = idxbx
        else:
            raise Exception('Invalid idxbx value. Exiting.')

    @lbu.setter
    def lbu(self, lbu):
        if type(lbu) == np.ndarray:
            self._lbu = lbu
        else:
            raise Exception('Invalid lbu value. Exiting.')

    @ubu.setter
    def ubu(self, ubu):
        if type(ubu) == np.ndarray:
            self._ubu = ubu
        else:
            raise Exception('Invalid ubu value. Exiting.')
    
    @idxbu.setter
    def idxbu(self, idxbu):
        if type(idxbu) == np.ndarray:
            self._idxbu = idxbu
        else:
            raise Exception('Invalid idxbu value. Exiting.')

    @lg.setter
    def lg(self, lg):
        if type(lg) == np.ndarray:
            self._lg = lg
        else:
            raise Exception('Invalid lg value. Exiting.')

    @ug.setter
    def ug(self, ug):
        if type(ug) == np.ndarray:
            self._ug = ug
        else:
            raise Exception('Invalid ug value. Exiting.')

    @lh.setter
    def lh(self, lh):
        if type(lh) == np.ndarray:
            self._lh = lh
        else:
            raise Exception('Invalid lh value. Exiting.')

    @uh.setter
    def uh(self, uh):
        if type(uh) == np.ndarray:
            self._uh = uh
        else:
            raise Exception('Invalid uh value. Exiting.')

    @D.setter
    def D(self, D):
        if type(D) == np.ndarray:
            self._D = D
        else:
            raise Exception('Invalid D value. Exiting.')

    @C.setter
    def C(self, C):
        if type(C) == np.ndarray:
            self._C = C
        else:
            raise Exception('Invalid C value. Exiting.')

    @CN.setter
    def CN(self, CN):
        if type(CN) == np.ndarray:
            self._CN = CN
        else:
            raise Exception('Invalid CN value. Exiting.')

    @lbxN.setter
    def lbxN(self, lbxN):
        if type(lbxN) == np.ndarray:
            self._lbxN = lbxN
        else:
            raise Exception('Invalid lbxN value. Exiting.')

    @ubxN.setter
    def ubxN(self, ubxN):
        if type(ubxN) == np.ndarray:
            self._ubxN = ubxN
        else:
            raise Exception('Invalid ubxN value. Exiting.')

    @idxbxN.setter
    def idxbxN(self, idxbxN):
        if type(idxbxN) == np.ndarray:
            self._idxbxN = idxbxN
        else:
            raise Exception('Invalid idxbxN value. Exiting.')

    @x0.setter
    def x0(self, x0):
        if type(x0) == np.ndarray:
            self._x0 = x0
        else:
            raise Exception('Invalid x0 value. Exiting.')

    @p.setter
    def p(self, p):
        if type(p) == np.ndarray:
            self._p = p
        else:
            raise Exception('Invalid p value. Exiting.')

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
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM')

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
        self.con_p_name = None 
        self.con_pN_name = None 
        self.con_h_name = None 
        self.con_hN_name = None 
        self.constants = []
        self.acados_include_path = []
        self.acados_lib_path = []

def check_ra(ra):
    if ra.solver_config.hessian_approx == 'EXACT' and ra.solver_config.integrator_type == 'IRK':
        raise Exception('Exact Hessians not yet supported with IRK integrators.')

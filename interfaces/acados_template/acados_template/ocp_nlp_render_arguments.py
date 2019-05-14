import numpy as np
import json
import os
import sys

class ocp_nlp_dims:
    """
    class containing the dimensions of the optimal control problem
    """
    def __init__(self):
        self.__nx    = None  #: number of states
        self.__nz    = 0     #: number of algebraic variables
        self.__nu    = None  #: number of inputs
        self.__np    = 0     #: number of parameters
        self.__ny    = None  #: number of residuals in Lagrange term
        self.__nyN   = None  #: number of residuals in Mayer term
        self.__npd   = 0     #: number of positive definite constraints
        self.__npdN  = 0     #: number of positive definite constraints at t=T
        self.__nh    = 0     #: number of nonlinear constraints
        self.__nhN   = 0     #: number of nonlinear constraints at t=T
        self.__nbx   = 0     #: number of state bounds 
        self.__nbxN  = 0     #: number of state bounds at t=T 
        self.__nbu   = 0     #: number of input bounds
        self.__nsbx  = 0     #: number of soft state bounds 
        self.__nsbxN = 0     #: number of soft state bounds at t=T 
        self.__nsbu  = 0     #: number of soft input bounds
        self.__ns    = 0     #: total number of slacks
        self.__nsN   = 0     #: total number of slacks at t=T
        self.__ng    = 0     #: number of general constraints
        self.__ngN   = 0     #: number of general constraints at t=T
        self.__N     = None  #: prediction horizon 

    @property
    def nx(self):
        return self.__nx

    @property
    def nz(self):
        return self.__nz

    @property
    def nu(self):
        return self.__nu

    @property
    def np(self):
        return self.__np

    @property
    def ny(self):
        return self.__ny

    @property
    def nyN(self):
        return self.__nyN

    @property
    def npd(self):
        return self.__npd

    @property
    def npdN(self):
        return self.__npdN

    @property
    def nh(self):
        return self.__nh

    @property
    def nhN(self):
        return self.__nhN

    @property
    def nbx(self):
        return self.__nbx

    @property
    def nbxN(self):
        return self.__nbxN

    @property
    def nbu(self):
        return self.__nbu

    @property
    def nsbx(self):
        return self.__nsbx

    @property
    def nsbxN(self):
        return self.__nsbx

    @property
    def nsbu(self):
        return self.__nsbu

    @property
    def ns(self):
        return self.__ns

    @property
    def nsN(self):
        return self.__nsN

    @property
    def ng(self):
        return self.__ng

    @property
    def ngN(self):
        return self.__ngN

    @property
    def N(self):
        return self.__N

    @nx.setter
    def nx(self, nx):
        if type(nx) == int and nx > 0:
            self.__nx = nx
        else:
            raise Exception('Invalid nx value. Exiting.')

    @nz.setter
    def nz(self, nz):
        if type(nz) == int and nz > -1:
            self.__nz = nz
        else:
            raise Exception('Invalid nz value. Exiting.')

    @nu.setter
    def nu(self, nu):
        if type(nu) == int and nu > 0:
            self.__nu = nu
        else:
            raise Exception('Invalid nu value. Exiting.')

    @np.setter
    def np(self, np):
        if type(np) == int and np > -1:
            self.__np = np
        else:
            raise Exception('Invalid np value. Exiting.')

    @ny.setter
    def ny(self, ny):
        if type(ny) == int and ny > -1:
            self.__ny = ny
        else:
            raise Exception('Invalid ny value. Exiting.')

    @nyN.setter
    def nyN(self, nyN):
        if type(nyN) == int and nyN > -1:
            self.__nyN = nyN
        else:
            raise Exception('Invalid nyN value. Exiting.')

    @npd.setter
    def npd(self, npd):
        if type(npd) == int and npd > -1:
            self.__npd = npd
        else:
            raise Exception('Invalid npd value. Exiting.')

    @npdN.setter
    def npdN(self, npdN):
        if type(npdN) == int and npdN > -1:
            self.__npdN = npdN
        else:
            raise Exception('Invalid npdN value. Exiting.')

    @nh.setter
    def nh(self, nh):
        if type(nh) == int and nh > -1:
            self.__nh = nh
        else:
            raise Exception('Invalid nh value. Exiting.')

    @nhN.setter
    def nhN(self, nhN):
        if type(nhN) == int and nhN > -1:
            self.__nhN = nhN
        else:
            raise Exception('Invalid nhN value. Exiting.')

    @nbx.setter
    def nbx(self, nbx):
        if type(nbx) == int and nbx > -1:
            self.__nbx = nbx
        else:
            raise Exception('Invalid nbx value. Exiting.')

    @nbxN.setter
    def nbxN(self, nbxN):
        if type(nbxN) == int and nbxN > -1:
            self.__nbxN = nbxN
        else:
            raise Exception('Invalid nbxN value. Exiting.')

    @nbu.setter
    def nbu(self, nbu):
        if type(nbu) == int and nbu > -1:
            self.__nbu = nbu
        else:
            raise Exception('Invalid nbu value. Exiting.')

    @nsbx.setter
    def nsbx(self, nbx):
        if type(nsbx) == int and nsbx > -1:
            self.__nsbx = nsbx
        else:
            raise Exception('Invalid nsbx value. Exiting.')

    @nsbxN.setter
    def nsbxN(self, nbxN):
        if type(nsbxN) == int and nsbxN > -1:
            self.__nsbxN = nsbxN
        else:
            raise Exception('Invalid nsbxN value. Exiting.')

    @nsbu.setter
    def nsbu(self, nsbu):
        if type(nsbu) == int and nsbu > -1:
            self.__nsbu = nsbu
        else:
            raise Exception('Invalid nsbu value. Exiting.')

    @ns.setter
    def ns(self, ns):
        if type(ns) == int and ns > -1:
            self.__ns = ns
        else:
            raise Exception('Invalid ns value. Exiting.')

    @nsN.setter
    def nsN(self, nsN):
        if type(nsN) == int and nsN > -1:
            self.__nsN = nsN
        else:
            raise Exception('Invalid nsN value. Exiting.')

    @ng.setter
    def ng(self, ng):
        if type(ng) == int and ng > -1:
            self.__ng = ng
        else:
            raise Exception('Invalid ng value. Exiting.')

    @ngN.setter
    def ngN(self, ngN):
        if type(ngN) == int and ngN > -1:
            self.__ngN = ngN
        else:
            raise Exception('Invalid ngN value. Exiting.')

    @N.setter
    def N(self, N):
        if type(N) == int and N > 0:
            self.__N = N
        else:
            raise Exception('Invalid N value. Exiting.')

class ocp_nlp_cost:
    """
    class containing the description of the cost
    """
    # linear least-squares cost: || Vx*x + Vu*x + Vz*z ||^2_W
    def __init__(self):
        # Lagrange term
        self.__W     = []  # weight matrix
        self.__Vx    = []  # x matrix coefficient
        self.__Vu    = []  # u matrix coefficient
        self.__Vz    = []  # z matrix coefficient
        self.__yref  = []  # reference
        self.__Zl    = []  # Hessian wrt lower slack 
        self.__Zu    = []  # Hessian wrt upper slack 
        self.__zl    = []  # gradient wrt lower slack 
        self.__zu    = []  # gradient wrt upper slack 
        # Mayer term
        self.__WN    = []  # weight matrix
        self.__VxN   = []  # x matrix coefficient
        self.__yrefN = []  # reference
        self.__ZlN   = []  # Hessian wrt lower slack 
        self.__ZuN   = []  # Hessian wrt upper slack 
        self.__zlN   = []  # gradient wrt lower slack 
        self.__zuN   = []  # gradient wrt upper slack 

    # Lagrange term
    @property
    def W(self):
        return self.__W

    @property
    def Vx(self):
        return self.__Vx

    @property
    def Vu(self):
        return self.__Vu

    @property
    def Vz(self):
        return self.__Vz

    @property
    def yref(self):
        return self.__yref

    @property
    def Zl(self):
        return self.__Zl

    @property
    def Zu(self):
        return self.__Zu

    @property
    def zl(self):
        return self.__zl

    @property
    def zu(self):
        return self.__zu

    @W.setter
    def W(self, W):
        if type(W) == np.ndarray:
            self.__W = W
        else:
            raise Exception('Invalid W value. Exiting.')
    
    @Vx.setter
    def Vx(self, Vx):
        if type(Vx) == np.ndarray:
            self.__Vx = Vx
        else:
            raise Exception('Invalid Vx value. Exiting.')
    
    @Vu.setter
    def Vu(self, Vu):
        if type(Vu) == np.ndarray:
            self.__Vu = Vu
        else:
            raise Exception('Invalid Vu value. Exiting.')

    @Vz.setter
    def Vz(self, Vz):
        if type(Vz) == np.ndarray:
            self.__Vz = Vz
        else:
            raise Exception('Invalid Vz value. Exiting.')

    @yref.setter
    def yref(self, yref):
        if type(yref) == np.ndarray:
            self.__yref = yref
        else:
            raise Exception('Invalid yref value. Exiting.')

    @Zl.setter
    def Zl(self, Zl):
        if type(Zl) == np.ndarray:
            self.__Zl = Zl
        else:
            raise Exception('Invalid Zl value. Exiting.')

    @Zu.setter
    def Zu(self, Zu):
        if type(Zu) == np.ndarray:
            self.__Zu = Zu
        else:
            raise Exception('Invalid Zu value. Exiting.')

    @zl.setter
    def zl(self, zl):
        if type(zl) == np.ndarray:
            self.__zl = zl
        else:
            raise Exception('Invalid zl value. Exiting.')

    @zu.setter
    def zu(self, zu):
        if type(zu) == np.ndarray:
            self.__zu = zu
        else:
            raise Exception('Invalid zu value. Exiting.')

    # Mayer term
    @property
    def WN(self):
        return self.__WN

    @property
    def VxN(self):
        return self.__VxN

    @property
    def yrefN(self):
        return self.__yrefN

    @property
    def ZlN(self):
        return self.__ZlN

    @property
    def ZuN(self):
        return self.__ZuN

    @property
    def zlN(self):
        return self.__zlN

    @property
    def zuN(self):
        return self.__zuN

    @WN.setter
    def WN(self, WN):
        if type(WN) == np.ndarray:
            self.__WN = WN
        else:
            raise Exception('Invalid WN value. Exiting.')
    
    @VxN.setter
    def VxN(self, VxN):
        if type(VxN) == np.ndarray:
            self.__VxN = VxN
        else:
            raise Exception('Invalid VxN value. Exiting.')

    @yrefN.setter
    def yrefN(self, yrefN):
        if type(yrefN) == np.ndarray:
            self.__yrefN = yrefN
        else:
            raise Exception('Invalid yrefN value. Exiting.')

    @ZlN.setter
    def ZlN(self, ZlN):
        if type(ZlN) == np.ndarray:
            self.__ZlN = ZlN
        else:
            raise Exception('Invalid ZlN value. Exiting.')

    @Zu.setter
    def ZuN(self, ZuN):
        if type(ZuN) == np.ndarray:
            self.__ZuN = ZuN
        else:
            raise Exception('Invalid ZuN value. Exiting.')

    @zlN.setter
    def zlN(self, zlN):
        if type(zlN) == np.ndarray:
            self.__zlN = zlN
        else:
            raise Exception('Invalid zlN value. Exiting.')

    @zuN.setter
    def zuN(self, zuN):
        if type(zuN) == np.ndarray:
            self.__zuN = zuN
        else:
            raise Exception('Invalid zuN value. Exiting.')

class ocp_nlp_constraints:
    """
    class containing the description of the constraints
    """
    def __init__(self):
        # bounds on x and u
        self.__lbx    = []  # lower bounds on x
        self.__lbu    = []  # lower bounds on u
        self.__ubx    = []  # upper bounds on x 
        self.__ubu    = []  # upper bounds on u 
        self.__idxbx  = []  # indexes of bounds on x 
        self.__idxbu  = []  # indexes of bounds on u
        # bounds on x at t=T
        self.__lbxN   = []  # lower bounds on x at t=T 
        self.__ubxN   = []  # upper bounds on x at t=T 
        self.__idxbxN = []  # indexes for bounds on x at t=T 
        # soft bounds on x and u
        self.__lsbx   = []  # soft lower bounds on x
        self.__lsbu   = []  # soft lower bounds on u
        self.__usbx   = []  # soft upper bounds on x 
        self.__usbu   = []  # soft upper bounds on u 
        self.__idxsbx = []  # indexes of soft bounds on x 
        self.__idxsbu = []  # indexes of soft bounds on u
        # soft bounds on x and u at t=T
        self.__lsbxN  = []  # soft lower bounds on x at t=T
        self.__usbxN  = []  # soft upper bounds on x at t=T
        self.__idxsbxN= []  # indexes of soft bounds on x at t=T 
        # polytopic constraints 
        self.__lg     = []  # lower bound for general inequalities 
        self.__ug     = []  # upper bound for general inequalities 
        self.__D      = []  # D matrix in lg <= D * u + C * x <= ug
        self.__C      = []  # C matrix in lg <= D * u + C * x <= ug
        # polytopic constraints at t=T 
        self.__CN     = []  # C matrix at t=T 
        self.__lgN    = []  # lower bound on general inequalities at t=T 
        self.__ugN    = []  # upper bound on general inequalities at t=T 
        # nonlinear constraints
        self.__lh     = []  # lower bound for nonlinear inequalities 
        self.__uh     = []  # upper bound for nonlinear inequalities 
        # nonlinear constraints at t=T
        self.__uhN    = []  # upper bound on nonlinear inequalities at t=T 
        self.__lhN    = []  # lower bound on nonlinear inequalities at t=T 
        self.__x0     = []  # initial state 
        self.__p      = []  # parameters 

    @property
    def lbx(self):
        return self.__lbx

    @property
    def lbu(self):
        return self.__lbu
    
    @property
    def ubx(self):
        return self.__ubx

    @property
    def ubu(self):
        return self.__ubu

    @property
    def idxbx(self):
        return self.__idxbx

    @property
    def idxbu(self):
        return self.__idxbu

    @property
    def lsbx(self):
        return self.__lsbx

    @property
    def lsbu(self):
        return self.__lsbu
    
    @property
    def usbx(self):
        return self.__usbx

    @property
    def usbu(self):
        return self.__usbu

    @property
    def idxsbx(self):
        return self.__idxsbx

    @property
    def idxsbu(self):
        return self.__idxsbu

    @property
    def lsbxN(self):
        return self.__lsbxN

    @property
    def usbxN(self):
        return self.__usbxN

    @property
    def idxsbxN(self):
        return self.__idxsbxN

    @property
    def lg(self):
        return self.__lg

    @property
    def ug(self):
        return self.__ug

    @property
    def lh(self):
        return self.__lh

    @property
    def uh(self):
        return self.__uh

    @property
    def D(self):
        return self.__D

    @property
    def C(self):
        return self.__C

    @property
    def lbxN(self):
        return self.__lbxN

    @property
    def ubxN(self):
        return self.__ubxN

    @property
    def idxbxN(self):
        return self.__idxbxN

    @property
    def CN(self):
        return self.__CN

    @property
    def lgN(self):
        return self.__lgN

    @property
    def ugN(self):
        return self.__ugN

    @property
    def lgN(self):
        return self.__lgN

    @property
    def ugN(self):
        return self.__ugN

    @property
    def x0(self):
        return self.__x0

    @property
    def p(self):
        return self.__p

    @lbx.setter
    def lbx(self, lbx):
        if type(lbx) == np.ndarray:
            self.__lbx = lbx
        else:
            raise Exception('Invalid lbx value. Exiting.')

    @ubx.setter
    def ubx(self, ubx):
        if type(ubx) == np.ndarray:
            self.__ubx = ubx
        else:
            raise Exception('Invalid ubx value. Exiting.')

    @idxbx.setter
    def idxbx(self, idxbx):
        if type(idxbx) == np.ndarray:
            self.__idxbx = idxbx
        else:
            raise Exception('Invalid idxbx value. Exiting.')

    @lbu.setter
    def lbu(self, lbu):
        if type(lbu) == np.ndarray:
            self.__lbu = lbu
        else:
            raise Exception('Invalid lbu value. Exiting.')

    @ubu.setter
    def ubu(self, ubu):
        if type(ubu) == np.ndarray:
            self.__ubu = ubu
        else:
            raise Exception('Invalid ubu value. Exiting.')
    
    @idxbu.setter
    def idxbu(self, idxbu):
        if type(idxbu) == np.ndarray:
            self.__idxbu = idxbu
        else:
            raise Exception('Invalid idxbu value. Exiting.')

    @lsbx.setter
    def lsbx(self, lsbx):
        if type(lsbx) == np.ndarray:
            self.__lsbx = lsbx
        else:
            raise Exception('Invalid lsbx value. Exiting.')

    @usbx.setter
    def usbx(self, usbx):
        if type(usbx) == np.ndarray:
            self.__usbx = usbx
        else:
            raise Exception('Invalid usbx value. Exiting.')

    @idxsbx.setter
    def idxsbx(self, idxsbx):
        if type(idxsbx) == np.ndarray:
            self.__idxsbx = idxsbx
        else:
            raise Exception('Invalid idxsbx value. Exiting.')

    @lsbu.setter
    def lsbu(self, lsbu):
        if type(lsbu) == np.ndarray:
            self.__lsbu = lsbu
        else:
            raise Exception('Invalid lsbu value. Exiting.')

    @usbu.setter
    def usbu(self, usbu):
        if type(usbu) == np.ndarray:
            self.__usbu = usbu
        else:
            raise Exception('Invalid usbu value. Exiting.')
    
    @idxsbu.setter
    def idxsbu(self, idxsbu):
        if type(idxsbu) == np.ndarray:
            self.__idxsbu = idxsbu
        else:
            raise Exception('Invalid idxsbu value. Exiting.')

    @lsbxN.setter
    def lsbxN(self, lsbxN):
        if type(lsbxN) == np.ndarray:
            self.__lsbxN = lsbxN
        else:
            raise Exception('Invalid lsbxN value. Exiting.')

    @usbxN.setter
    def usbxN(self, usbxN):
        if type(usbxN) == np.ndarray:
            self.__usbxN = usbxN
        else:
            raise Exception('Invalid usbxN value. Exiting.')

    @idxsbxN.setter
    def idxsbxN(self, idxsbxN):
        if type(idxsbxN) == np.ndarray:
            self.__idxsbxN = idxsbxN
        else:
            raise Exception('Invalid idxsbxN value. Exiting.')

    @lg.setter
    def lg(self, lg):
        if type(lg) == np.ndarray:
            self.__lg = lg
        else:
            raise Exception('Invalid lg value. Exiting.')

    @ug.setter
    def ug(self, ug):
        if type(ug) == np.ndarray:
            self.__ug = ug
        else:
            raise Exception('Invalid ug value. Exiting.')

    @lh.setter
    def lh(self, lh):
        if type(lh) == np.ndarray:
            self.__lh = lh
        else:
            raise Exception('Invalid lh value. Exiting.')

    @uh.setter
    def uh(self, uh):
        if type(uh) == np.ndarray:
            self.__uh = uh
        else:
            raise Exception('Invalid uh value. Exiting.')

    @D.setter
    def D(self, D):
        if type(D) == np.ndarray:
            self.__D = D
        else:
            raise Exception('Invalid D value. Exiting.')

    @C.setter
    def C(self, C):
        if type(C) == np.ndarray:
            self.__C = C
        else:
            raise Exception('Invalid C value. Exiting.')

    @CN.setter
    def CN(self, CN):
        if type(CN) == np.ndarray:
            self.__CN = CN
        else:
            raise Exception('Invalid CN value. Exiting.')

    @lbxN.setter
    def lbxN(self, lbxN):
        if type(lbxN) == np.ndarray:
            self.__lbxN = lbxN
        else:
            raise Exception('Invalid lbxN value. Exiting.')

    @ubxN.setter
    def ubxN(self, ubxN):
        if type(ubxN) == np.ndarray:
            self.__ubxN = ubxN
        else:
            raise Exception('Invalid ubxN value. Exiting.')

    @idxbxN.setter
    def idxbxN(self, idxbxN):
        if type(idxbxN) == np.ndarray:
            self.__idxbxN = idxbxN
        else:
            raise Exception('Invalid idxbxN value. Exiting.')

    @x0.setter
    def x0(self, x0):
        if type(x0) == np.ndarray:
            self.__x0 = x0
        else:
            raise Exception('Invalid x0 value. Exiting.')

    @p.setter
    def p(self, p):
        if type(p) == np.ndarray:
            self.__p = p
        else:
            raise Exception('Invalid p value. Exiting.')

class ocp_nlp_solver_config:
    """
    class containing the description of the solver configuration
    """
    def __init__(self):
        self.__qp_solver        = 'PARTIAL_CONDENSING_HPIPM'  # qp solver to be used in the NLP solver
        self.__hessian_approx   = 'GAUSS_NEWTON'              # hessian approximation
        self.__integrator_type  = 'ERK'                       # integrator type
        self.__tf               = None                        # prediction horizon
        self.__nlp_solver_type  = 'SQP_RTI'                   # NLP solver 

    @property
    def qp_solver(self):
        return self.__qp_solver

    @property
    def hessian_approx(self):
        return self.__hessian_approx

    @property
    def integrator_type(self):
        return self.__integrator_type

    @property
    def nlp_solver_type(self):
        return self.__nlp_solver_type

    @qp_solver.setter
    def qp_solver(self, qp_solver):
        qp_solvers = ('PARTIAL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_QPOASES', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM')

        if type(qp_solver) == str and qp_solver in qp_solvers:
            self.__qp_solver = qp_solver
        else:
            raise Exception('Invalid qp_solver value. Possible values are:\n\n' \
                    + ',\n'.join(qp_solvers) + '.\n\nYou have: ' + qp_solver + '.\n\nExiting.')
    @property
    def tf(self):
        return self.__tf

    @hessian_approx.setter
    def hessian_approx(self, hessian_approx):
        hessian_approxs = ('GAUSS_NEWTON')

        if type(hessian_approx) == str and hessian_approx in hessian_approxs:
            self.__hessian_approx = hessian_approx
        else:
            raise Exception('Invalid hessian_approx value. Possible values are:\n\n' \
                    + ',\n'.join(hessian_approxs) + '.\n\nYou have: ' + hessian_approx + '.\n\nExiting.')

    @integrator_type.setter
    def integrator_type(self, integrator_type):
        integrator_types = ('ERK', 'IRK')

        if type(integrator_type) == str and integrator_type in integrator_types:
            self.__integrator_type = integrator_type
        else:
            raise Exception('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(integrator_types) + '.\n\nYou have: ' + integrator_type + '.\n\nExiting.')

    @tf.setter
    def tf(self, tf):
        self.__tf = tf

    @nlp_solver_type.setter
    def nlp_solver_type(self, nlp_solver_type):
        nlp_solver_types = ('SQP', 'SQP_RTI')

        if type(nlp_solver_type) == str and nlp_solver_type in nlp_solver_types:
            self.__nlp_solver_type = nlp_solver_type
        else:
            raise Exception('Invalid nlp_solver_type value. Possible values are:\n\n' \
                    + ',\n'.join(nlp_solver_types) + '.\n\nYou have: ' + nlp_solver_type + '.\n\nExiting.')

class acados_ocp_nlp:
    """
    class containing the full description if the optimal control problem
    """
    def __init__(self):
        self.dims = ocp_nlp_dims()
        self.cost = ocp_nlp_cost()
        self.constraints = ocp_nlp_constraints()
        self.solver_config = ocp_nlp_solver_config()
        self.model_name  = None 
        self.con_p_name  = None 
        self.con_pN_name = None 
        self.con_h_name  = None 
        self.con_hN_name = None 
        self.constants = {}
        self.acados_include_path = []
        self.acados_lib_path = []

def check_ra(ra):
    """
    (DEPRECATED) function that checks the consistency of the optimal control description
    """
    # TODO(andrea): dimensions check are already performed 
    # on the JSON data and type checks should be enforced by the 
    # property setters. Add extra checks here?
    return

def np_array_to_list(np_array):
    return np_array.tolist()

class ocp_nlp_as_object:
        def __init__(self, d):
            self.__dict__ = d

def dict2json(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict2json(v)

        v_type = str(type(v).__name__)
        # out_key = '__' + v_type + '__' + k.split('__', 1)[-1]
        out_key = k.split('__', 1)[-1]
        out[k.replace(k, out_key)] = v
    return out

def acados_ocp2json_layout(acados_ocp):
    """ Convert acados ocp nlp object JSON format by stripping the 
    property mangling and adding array dimension info.
    ALL items of type String will be converted 
    to type ndarrray!
     
    Parameters
    ----------
    acados_ocp : class
        object of type ocp_nlp_render_arguments.
    
    Returns
    ------
    out: dict 
        acados_layout
    """
    ocp_nlp = acados_ocp
    ocp_nlp.cost = acados_ocp.cost.__dict__
    ocp_nlp.constraints = acados_ocp.constraints.__dict__
    ocp_nlp.solver_config = acados_ocp.solver_config.__dict__
    ocp_nlp.dims = acados_ocp.dims.__dict__
    ocp_nlp = ocp_nlp.__dict__
    json_layout = dict2json_layout(ocp_nlp)
    return json_layout

def dict2json_layout(d):
    """ Convert dictionary containing the description of 
    of the ocp_nlp to JSON format by stripping the 
    property mangling and adding array dimension info.
    ALL items of type String will be converted 
    to type ndarrray!
     
    Parameters
    ----------
    d : dict
        dictionary containing the description of 
        the ocp_nlp.
    
    Returns
    ------
    out: dict 
        postprocessed dictionary.
    """
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict2json_layout(v)

        v_type = str(type(v).__name__)
        if v_type == 'list':
            v_type = 'ndarray'

        # add array number of dimensions?
        # if v_type == 'ndarray':
        #     v_type = v_type + '_' + str(len(v.shape))
        out_key = k.split('__', 1)[-1]

        if isinstance(v, dict):
            out[k.replace(k, out_key)] = v  
        else:
            out[k.replace(k, out_key)] = [v_type] 
    
    return out

def cast_ocp_nlp(ocp_nlp, ocp_nlp_layout):
    """ MATLAB does not allow distinction between e.g a = [1,1,1] and b = [1,1,1].' 
    or a = 1 and b = [1]. Hence, we need to do some postprocessing of the JSON 
    file generated from MATLAB.
     
    Parameters
    ----------
    ocp_nlp : dict
        ocp_nlp dictionary to be postprocessed.
    
    ocp_nlp_layout : dict
        acados ocp_nlp target layout
    Returns
    ------
    out : dict
        postprocessed dictionary
    """

    out = {}
    for k, v in ocp_nlp.items():
        if isinstance(v, dict):
            v = cast_ocp_nlp(v, ocp_nlp_layout[k])

        if 'ndarray' in ocp_nlp_layout[k]:
            if isinstance(v, int) or isinstance(v, float):
                v = np.array([v])
        out[k] = v
    return out 

def json2dict(ocp_nlp, ocp_nlp_dims):
    # load JSON layout
    current_module = sys.modules[__name__]
    acados_path = os.path.dirname(current_module.__file__)
    with open(acados_path + '/acados_layout.json', 'r') as f:
        ocp_nlp_layout = json.load(f)

    out = json2dict_rec(ocp_nlp, ocp_nlp_dims, ocp_nlp_layout)
    return out

def json2dict_rec(ocp_nlp, ocp_nlp_dims, ocp_nlp_layout):
    """ convert ocp_nlp loaded JSON to dictionary. Mainly convert
    lists to arrays for easier handling.
    Parameters
    ---------
    ocp_nlp : dict 
        dictionary loaded from JSON to be post-processed.
    
    ocp_nlp_dims : dict 
        dictionary containing the ocp_nlp dimensions.

    ocp_nlp_layout : dict 
        acados ocp_nlp layout.

    Returns
    -------
    out : dict 
        post-processed dictionary.
    """
    out = {}
    for k, v in ocp_nlp.items():
        if isinstance(v, dict):
            v = json2dict_rec(v, ocp_nlp_dims, ocp_nlp_layout[k])

        v_type__ = str(type(v).__name__)
        out_key = k.split('__', 1)[-1]
        v_type = out_key.split('__')[0]
        out_key = out_key.split('__', 1)[-1]
        if 'ndarray' in ocp_nlp_layout[k]:
            if isinstance(v, int) or isinstance(v, float):
                v = np.array([v])
        if v_type == 'ndarray' or v_type__ == 'list':
            if v == []:
                # v = None
                v = []
            else:
                v = np.array(v)
                v_dims = v.shape
                dim_keys = ocp_nlp_layout[k][1]
                dims_l = []
                dims_names = []
                for item in dim_keys:
                    dims_l.append(ocp_nlp_dims[item])
                    dims_names.append(item)
                dims = tuple(dims_l)
                try: 
                    v = np.reshape(v, dims)
                except:  
                    raise Exception('acados -- mismatching dimensions for field {0}. Provided data has dimensions {1}, while associated dimensions {2} are {3}'.format(out_key, v_dims, dims_names, dims))
        out[k.replace(k, out_key)] = v
    return out

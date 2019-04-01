import numpy as np

class ocp_nlp_dims:
    def __init__(self):
        self.__nx   = None  # number of states
        self.__nz   = 0     # number of algebraic variables
        self.__nu   = None  # number of inputs
        self.__np   = 0     # number of parameters
        self.__ny   = None  # number of residuals in Lagrange term
        self.__nyN  = None  # number of residuals in Mayer term
        self.__npd  = 0     # number of positive definite constraints
        self.__npdN = 0     # number of positive definite constraints in last stage
        self.__nh   = 0     # number of nonlinear constraints
        self.__nhN  = 0     # number of nonlinear constraints in last stage
        self.__nbx  = 0     # number of state bounds 
        self.__nbu  = 0     # number of input bounds
        self.__ng   = 0     # number of general constraints
        self.__nbxN = 0     # number of state bounds in last stage 
        self.__ngN  = 0     # number of general constraints in last stage
        self.__N    = None  # prediction horizon 

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
    def nyN(self):
        return self.__nyN

    @property
    def nbx(self):
        return self.__nbx

    @property
    def nbu(self):
        return self.__nbu

    @property
    def ng(self):
        return self.__ng

    @property
    def nbxN(self):
        return self.__nbxN
    
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

    @nbu.setter
    def nbu(self, nbu):
        if type(nbu) == int and nbu > -1:
            self.__nbu = nbu
        else:
            raise Exception('Invalid nbu value. Exiting.')

    @nbx.setter
    def nbx(self, nbx):
        if type(nbx) == int and nbx > -1:
            self.__nbx = nbx
        else:
            raise Exception('Invalid nbx value. Exiting.')

    @ng.setter
    def ng(self, ng):
        if type(ng) == int and ng > -1:
            self.__ng = ng
        else:
            raise Exception('Invalid ng value. Exiting.')

    @nbxN.setter
    def nbxN(self, nbxN):
        if type(nbxN) == int and nbxN > -1:
            self.__nbxN = nbxN
        else:
            raise Exception('Invalid nbxN value. Exiting.')

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
    # linear least-squares cost: || Vx*x + Vu*x + Vz*z ||^2_W
    def __init__(self):
        # Lagrange term
        self.__W     = None  # weight matrix
        self.__Vx    = None  # x matrix coefficient
        self.__Vu    = None  # u matrix coefficient
        self.__Vz    = None  # z matrix coefficient
        self.__yref  = None  # reference
        # Mayer term
        self.__WN    = None  # weight matrix
        self.__VxN   = None  # x matrix coefficient
        self.__yrefN = None  # reference

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

class ocp_nlp_constraints:
    def __init__(self):
        self.__lbx    = None  # lower bounds on x
        self.__lbu    = None  # lower bounds on u
        self.__idxbx  = None  # indexes of bounds on x 
        self.__ubx    = None  # upper bounds on x 
        self.__ubu    = None  # upper bounds on u 
        self.__idxbu  = None  # indexes of bounds on u
        self.__lg     = None  # lower bound for general inequalities 
        self.__ug     = None  # upper bound for general inequalities 
        self.__lh     = None  # lower bound for nonlinear inequalities 
        self.__uh     = None  # upper bound for nonlinear inequalities 
        self.__D      = None  # D matrix in lg <= D * u + C * x <= ug
        self.__C      = None  # C matrix in lg <= D * u + C * x <= ug
        self.__lbxN   = None  # lower bounds on x at t=T 
        self.__ubxN   = None  # upper bounds on x at t=T 
        self.__idxbxN = None  # indexes for bounds on x at t=T 
        self.__CN     = None  # C matrix at t=T 
        self.__lgN    = None  # lower bound on general inequalities at t=T 
        self.__ugN    = None  # upper bound on general inequalities at t=T 
        self.__lhN    = None  # lower bound on nonlinear inequalities at t=T 
        self.__uhN    = None  # upper bound on nonlinear inequalities at t=T 
        self.__x0     = None  # initial state 

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
        hessian_approxs = ('GAUSS_NEWTON', 'EXACT')

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
        self.constants = {}
        self.acados_include_path = []
        self.acados_lib_path = []

def check_ra(ra):
    if ra.solver_config.hessian_approx == 'EXACT' and ra.solver_config.integrator_type == 'IRK':
        raise Exception('Exact Hessians not yet supported with IRK integrators.')

def np_array_to_list(np_array):
    return np_array.tolist()

class ocp_nlp_as_object:
        def __init__(self, d):
            self.__dict__ = d

def dict2json(d):
    new = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict2json(v)

        v_type = str(type(v).__name__)
        new_key = '__' + v_type + '__' + k.split('__', 1)[-1]
        new[k.replace(k, new_key)] = v
    return new

def dict2json_layout(d):
    new = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict2json_layout(v)

        v_type = str(type(v).__name__)
        if v_type == 'ndarray':
            v_type = v_type + '_' + str(len(v.shape))
        new_key = k.split('__', 1)[-1]

        if isinstance(v, dict):
            new[k.replace(k, new_key)] = v  
        else:
            new[k.replace(k, new_key)] = v_type 
    
    return new

def cast_ocp_nlp(ocp_nlp, layout):

    
    return ocp_nlp

def cast_ocp_nlp(ocp_nlp, ocp_nlp_layout):
    new = {}
    for k, v in ocp_nlp.items():
        if isinstance(v, dict):
            v = cast_ocp_nlp(v, ocp_nlp_layout[k])

        if 'ndarray' in ocp_nlp_layout[k]:
            dim = int(ocp_nlp_layout[k].split('_',1)[1])
            if isinstance(v, int) or isinstance(v, float):
                import pdb; pdb.set_trace()
                v = [v]
            print(v)
            while len(v.shape) < dim: 
                v = [v]
        new[k] = v
    return new 

    # load JSON layout
    with open('acados_layout.json', 'r') as f:
        ocp_nlp_layout = json.load(f)
    
    return ocp_nlp

def json2dict(d):
    new = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = json2dict(v)

        v_type__ = str(type(v).__name__)
        new_key = k.split('__', 1)[-1]
        v_type = new_key.split('__')[0]
        new_key = new_key.split('__', 1)[-1]
        # TODO: cast v to corresponding type
        print('key = ', k)
        print('v type = ', v_type__)
        print('v = ', v, '\n')
        if v_type == 'ndarray' or v_type__ == 'list':
            if v == []:
                v = None
            else:
                v = np.array(v)

        new[k.replace(k, new_key)] = v
    return new

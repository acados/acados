import numpy as np
import json
import os
import sys
from .casadi_functions import *

class ocp_nlp_dims:
    """
    class containing the dimensions of the optimal control problem
    """
    def __init__(self):
        self.__nx     = None  #: :math:`n_x` - number of states 
        self.__nz     = 0     #: :math:`n_z` - number of algebraic variables 
        self.__nu     = None  #: :math:`n_u` - number of inputs 
        self.__np     = 0     #: :math:`n_p` - number of parameters 
        self.__ny     = None  #: :math:`n_y` - number of residuals in Lagrange term 
        self.__ny_e   = None  #: :math:`n_{y}^e` - number of residuals in Mayer term 
        self.__npd    = 0     #: :math:`n_{\pi}` - dimension of the image of the inner nonlinear function in positive definite constraints 
        self.__npd_e  = 0     #: :math:`n_{\pi}^e` - dimension of the image of the inner nonlinear function in positive definite constraints
        self.__nh     = 0     #: :math:`n_h` - number of nonlinear constraints 
        self.__nh_e   = 0     #: :math:`n_{h}^e` - number of nonlinear constraints at t=T 
        self.__nbx    = 0     #: :math:`n_{b_x}` - number of state bounds 
        self.__nbx_e  = 0     #: :math:`n_{b_x}` - number of state bounds at t=T 
        self.__nbu    = 0     #: :math:`n_{b_u}` - number of input bounds 
        self.__nsbx   = 0     #: :math:`n_{{sb}_x}` - number of soft state bounds 
        self.__nsbx_e = 0     #: :math:`n_{{sb}^e_{x}}` - number of soft state bounds at t=T 
        self.__nsbu   = 0     #: :math:`n_{{sb}_u}` - number of soft input bounds 
        self.__nsh    = 0     #: :math:`n_{{sb}_u}` - number of soft nonlinear constraints 
        self.__nsh_e  = 0     #: :math:`n_{{sb}_u}` - number of soft nonlinear constraints 
        self.__ns     = 0     #: :math:`n_{s}` - total number of slacks 
        self.__ns_e   = 0     #: :math:`n_{s}^e` - total number of slacks at t=T 
        self.__ng     = 0     #: :math:`n_{g}` - number of general polytopic constraints 
        self.__ng_e   = 0     #: :math:`n_{g}^e` - number of general polytopic constraints at t=T 
        self.__N      = None  #: :math:`N` - prediction horizon  

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
    def ny_e(self):
        return self.__ny_e

    @property
    def npd(self):
        return self.__npd

    @property
    def npd_e(self):
        return self.__npd_e

    @property
    def nh(self):
        return self.__nh

    @property
    def nh_e(self):
        return self.__nh_e

    @property
    def nbx(self):
        return self.__nbx

    @property
    def nbx_e(self):
        return self.__nbx_e

    @property
    def nbu(self):
        return self.__nbu

    @property
    def nsbx(self):
        return self.__nsbx

    @property
    def nsbx_e(self):
        return self.__nsbx

    @property
    def nsbu(self):
        return self.__nsbu

    @property
    def nsh(self):
        return self.__nsh

    @property
    def nsh_e(self):
        return self.__nsh_e

    @property
    def ns(self):
        return self.__ns

    @property
    def ns_e(self):
        return self.__ns_e

    @property
    def ng(self):
        return self.__ng

    @property
    def ng_e(self):
        return self.__ng_e

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

    @ny_e.setter
    def ny_e(self, ny_e):
        if type(ny_e) == int and ny_e > -1:
            self.__ny_e = ny_e
        else:
            raise Exception('Invalid ny_e value. Exiting.')

    @npd.setter
    def npd(self, npd):
        if type(npd) == int and npd > -1:
            self.__npd = npd
        else:
            raise Exception('Invalid npd value. Exiting.')

    @npd_e.setter
    def npd_e(self, npd_e):
        if type(npd_e) == int and npd_e > -1:
            self.__npd_e = npd_e
        else:
            raise Exception('Invalid npd_e value. Exiting.')

    @nh.setter
    def nh(self, nh):
        if type(nh) == int and nh > -1:
            self.__nh = nh
        else:
            raise Exception('Invalid nh value. Exiting.')

    @nh_e.setter
    def nh_e(self, nh_e):
        if type(nh_e) == int and nh_e > -1:
            self.__nh_e = nh_e
        else:
            raise Exception('Invalid nh_e value. Exiting.')

    @nbx.setter
    def nbx(self, nbx):
        if type(nbx) == int and nbx > -1:
            self.__nbx = nbx
        else:
            raise Exception('Invalid nbx value. Exiting.')

    @nbx_e.setter
    def nbx_e(self, nbx_e):
        if type(nbx_e) == int and nbx_e > -1:
            self.__nbx_e = nbx_e
        else:
            raise Exception('Invalid nbx_e value. Exiting.')

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

    @nsbx_e.setter
    def nsbx_e(self, nbx_e):
        if type(nsbx_e) == int and nsbx_e > -1:
            self.__nsbx_e = nsbx_e
        else:
            raise Exception('Invalid nsbx_e value. Exiting.')

    @nsbu.setter
    def nsbu(self, nsbu):
        if type(nsbu) == int and nsbu > -1:
            self.__nsbu = nsbu
        else:
            raise Exception('Invalid nsbu value. Exiting.')

    @nsh.setter
    def nsh(self, nsh):
        if type(nsh) == int and nsh > -1:
            self.__nsh = nsh
        else:
            raise Exception('Invalid nsh value. Exiting.')

    @nsh_e.setter
    def nsh_e(self, nsh_e):
        if type(nsh_e) == int and nsh_e > -1:
            self.__nsh_e = nsh_e
        else:
            raise Exception('Invalid nsh_e value. Exiting.')

    @ns.setter
    def ns(self, ns):
        if type(ns) == int and ns > -1:
            self.__ns = ns
        else:
            raise Exception('Invalid ns value. Exiting.')

    @ns_e.setter
    def ns_e(self, ns_e):
        if type(ns_e) == int and ns_e > -1:
            self.__ns_e = ns_e
        else:
            raise Exception('Invalid ns_e value. Exiting.')

    @ng.setter
    def ng(self, ng):
        if type(ng) == int and ng > -1:
            self.__ng = ng
        else:
            raise Exception('Invalid ng value. Exiting.')

    @ng_e.setter
    def ng_e(self, ng_e):
        if type(ng_e) == int and ng_e > -1:
            self.__ng_e = ng_e
        else:
            raise Exception('Invalid ng_e value. Exiting.')

    @N.setter
    def N(self, N):
        if type(N) == int and N > 0:
            self.__N = N
        else:
            raise Exception('Invalid N value. Exiting.')

    def set(self, attr, value):
        setattr(self, attr, value)

class ocp_nlp_cost:
    """
    class containing the description of the cost
    (linear least-squares cost for the time being) 
    :math:`l(x,u,z) = || V_x x + V_u u + V_z z - y_{\\text{ref}}||^2_W`, 
    :math:`m(x) = || V^e_x x - y_{\\text{ref}^e}||^2_{W^e}`
    """
    def __init__(self):
        # Lagrange term
        self.__cost_type   = 'LINEAR_LS'  #: cost type
        self.__W           = []           #: :math:`W` - weight matrix
        self.__Vx          = []           #: :math:`V_x` - x matrix coefficient
        self.__Vu          = []           #: :math:`V_u` - u matrix coefficient
        self.__Vz          = []           #: :math:`V_z` - z matrix coefficient
        self.__yref        = []           #: :math:`y_{\text{ref}}` - reference
        self.__Zl          = []           #: :math:`Z_l` - Hessian wrt lower slack 
        self.__Zu          = []           #: :math:`Z_u` - Hessian wrt upper slack 
        self.__zl          = []           #: :math:`z_l` - gradient wrt lower slack 
        self.__zu          = []           #: :math:`z_u` - gradient wrt upper slack 
        # Mayer term
        self.__cost_type_e = 'LINEAR_LS'  #: cost type for Mayer term
        self.__W_e         = []           #: :math:`W^e` - weight matrix for Mayer term
        self.__Vx_e        = []           #: :math:`V_x^e` - x matrix coefficient for Mayer term
        self.__yref_e      = []           #: :math:`y_{\text{ref}}^e` - reference for Mayer term
        self.__Zl_e        = []           #: :math:`Z_l^e` - Hessian wrt lower slack for Mayer term
        self.__Zu_e        = []           #: :math:`Z_u^e` - Hessian wrt upper slack for Mayer term
        self.__zl_e        = []           #: :math:`z_l^e` - gradient wrt lower slack for Mayer term
        self.__zu_e        = []           #: :math:`z_u^e` - gradient wrt upper slack for Mayer term

    # Lagrange term
    @property
    def cost_type(self):
        return self.__cost_type

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

    @cost_type.setter
    def cost_type(self, cost_type):
        cost_types = ('LINEAR_LS')

        if type(cost_type) == str and cost_type in cost_types:
            self.__cost_type = cost_type
        else:
            raise Exception('Invalid cost_type value. Exiting.')

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
    def cost_type_e(self):
        return self.__cost_type_e

    @property
    def W_e(self):
        return self.__W_e

    @property
    def Vx_e(self):
        return self.__Vx_e

    @property
    def yref_e(self):
        return self.__yref_e

    @property
    def Zl_e(self):
        return self.__Zl_e

    @property
    def Zu_e(self):
        return self.__Zu_e

    @property
    def zl_e(self):
        return self.__zl_e

    @property
    def zu_e(self):
        return self.__zu_e

    @cost_type_e.setter
    def cost_type_e(self, cost_type_e):
        cost_types = ('LINEAR_LS')

        if type(cost_type_e) == str and cost_type_e in cost_types:
            self.__cost_type_e = cost_type_e
        else:
            raise Exception('Invalid cost_type_e value. Exiting.')

    @W_e.setter
    def W_e(self, W_e):
        if type(W_e) == np.ndarray:
            self.__W_e = W_e
        else:
            raise Exception('Invalid W_e value. Exiting.')
    
    @Vx_e.setter
    def Vx_e(self, Vx_e):
        if type(Vx_e) == np.ndarray:
            self.__Vx_e = Vx_e
        else:
            raise Exception('Invalid Vx_e value. Exiting.')

    @yref_e.setter
    def yref_e(self, yref_e):
        if type(yref_e) == np.ndarray:
            self.__yref_e = yref_e
        else:
            raise Exception('Invalid yref_e value. Exiting.')

    @Zl_e.setter
    def Zl_e(self, Zl_e):
        if type(Zl_e) == np.ndarray:
            self.__Zl_e = Zl_e
        else:
            raise Exception('Invalid Zl_e value. Exiting.')

    @Zu_e.setter
    def Zu_e(self, Zu_e):
        if type(Zu_e) == np.ndarray:
            self.__Zu_e = Zu_e
        else:
            raise Exception('Invalid Zu_e value. Exiting.')

    @zl_e.setter
    def zl_e(self, zl_e):
        if type(zl_e) == np.ndarray:
            self.__zl_e = zl_e
        else:
            raise Exception('Invalid zl_e value. Exiting.')

    @zu_e.setter
    def zu_e(self, zu_e):
        if type(zu_e) == np.ndarray:
            self.__zu_e = zu_e
        else:
            raise Exception('Invalid zu_e value. Exiting.')

    def set(self, attr, value):
        setattr(self, attr, value)

class ocp_nlp_constraints:
    """
    class containing the description of the constraints
    """
    def __init__(self):
        self.__constr_type   = 'BGH' #: constraint type
        self.__constr_type_e = 'BGH' #: constraint type
        # bounds on x and u
        self.__lbx     = []        #: :math:`\underline{x}` - lower bounds on x
        self.__lbu     = []        #: :math:`\underline{u}` - lower bounds on u
        self.__ubx     = []        #: :math:`\bar{x}` - upper bounds on x 
        self.__ubu     = []        #: :math:`\bar{u}` - upper bounds on u 
        self.__idxbx   = []        #: indexes of bounds on x (defines :math:`\Pi_x`) 
        self.__Jbx     = []        #: :math`J_x` - matrix coefficient for bounds on x 
        self.__idxbu   = []        #: indexes of bounds on u (defines :math:`\Pi_u`)
        self.__Jbu     = []        #: :math`J_u` - matrix coefficient for bounds on u 
        # bounds on x at t=T
        self.__lbx_e   = []        #: :math:`\underline{x}^e` - lower bounds on x at t=T 
        self.__ubx_e   = []        #: :math:`\bar{x}^e` - upper bounds on x at t=T 
        self.__idxbx_e = []        #: indexes for bounds on x at t=T (defines :math:`\Pi_x^e`) 
        self.__Jbx_e   = []         #: :math`J_{x}^e`indexes of bounds on x (defines :math:`\Pi_x`) 
        # polytopic constraints 
        self.__lg      = []        #: :math:`\underline{c}` - lower bound for general polytopic inequalities 
        self.__ug      = []        #: :math:`\bar{c}` - upper bound for general polytopic inequalities 
        self.__D       = []        #: :math:`D` - D matrix in lg <= D * u + C * x <= ug
        self.__C       = []        #: :math:`C` - C matrix in lg <= D * u + C * x <= ug
        # polytopic constraints at t=T 
        self.__C_e     = []        #: :math:`C^e` - C matrix at t=T 
        self.__lg_e    = []        #: :math:`\underline{c}^e` - lower bound on general polytopic inequalities at t=T 
        self.__ug_e    = []        #: :math:`\bar{c}^e` - upper bound on general polytopic inequalities at t=T 
        # nonlinear constraints
        self.__lh      = []        #: :math:`\underline{h}` - lower bound for nonlinear inequalities 
        self.__uh      = []        #: :math:`\bar{h}` - upper bound for nonlinear inequalities 
        # nonlinear constraints at t=T
        self.__uh_e    = []        #: :math:`\bar{h}^e` - upper bound on nonlinear inequalities at t=T 
        self.__lh_e    = []        #: :math:`\underline{h}^e` - lower bound on nonlinear inequalities at t=T 
        # soft bounds on x and u
        self.__lsbx   = []         #: soft lower bounds on x
        self.__lsbu   = []         #: soft lower bounds on u
        self.__usbx   = []         #: soft upper bounds on x 
        self.__usbu   = []         #: soft upper bounds on u 
        self.__idxsbx = []         #: indexes of soft bounds on x 
        self.__Jsbx   = []         #: :math`J_{s,x}` - matrix coefficient for soft bounds on x 
        self.__idxsbu = []         #: indexes of soft bounds on u
        self.__Jsbu   = []         #: :math`J_{s,u}` - matrix coefficient for soft bounds on u 
        # soft bounds on x at t=T
        self.__lsbx_e  = []        #: soft lower bounds on x at t=T
        self.__usbx_e  = []        #: soft upper bounds on x at t=T
        self.__idxsbx_e= []        #: indexes of soft bounds on x at t=T 
        self.__Jsbx_e    = []      #: :math`J_{s,x}^e` - matrix coefficient for soft bounds on x at t=T 
        # soft bounds on nonlinear constraints
        self.__lsh    = []         #: soft lower bounds for nonlinear constraints 
        self.__ush    = []         #: soft upper bounds for nonlinear constraints 
        self.__idxsh  = []         #: indexes of soft nonlinear constraints 
        self.__Jsh    = []         #: :math`J_{s,h}` - matrix coefficient for soft bounds on nonlinear constraints 
        # soft bounds on nonlinear constraints at t=T
        self.__lsh_e    = []       #: soft lower bounds for nonlinear constraints at t=T
        self.__ush_e    = []       #: soft upper bounds for nonlinear constraints at t=T
        self.__idxsh_e  = []       #: indexes of soft nonlinear constraints at t=T 
        self.__Jsh_e    = []       #: :math`J_{s,h}^e` - matrix coefficient for soft bounds on nonlinear constraints at t=T 
        self.__x0      = []        #: :math:`\bar{x}_0` - initial state 
        self.__p       = []        #: :math:`p` - parameters 

    # bounds on x and u
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
    def Jbx(self):
        return self.__Jbx

    @property
    def idxbu(self):
        return self.__idxbu

    @property
    def Jbu(self):
        return self.__Jbu

    # bounds on x at t=T
    @property
    def lbx_e(self):
        return self.__lbx_e

    @property
    def ubx_e(self):
        return self.__ubx_e

    @property
    def idxbx_e(self):
        return self.__idxbx_e

    @property
    def Jbx_e(self):
        return self.__Jbx_e

    # polytopic constraints 
    @property
    def C(self):
        return self.__C

    @property
    def D(self):
        return self.__D

    @property
    def lg(self):
        return self.__lg

    @property
    def ug(self):
        return self.__ug

    # polytopic constraints at t=T 
    @property
    def C_e(self):
        return self.__C_e

    @property
    def lg_e(self):
        return self.__lg_e

    @property
    def ug_e(self):
        return self.__ug_e

    # nonlinear constraints
    @property
    def lh(self):
        return self.__lh

    @property
    def uh(self):
        return self.__uh

    # nonlinear constraints at t=T
    @property
    def lh_e(self):
        return self.__lh_e

    @property
    def uh_e(self):
        return self.__uh_e

    # soft bounds on x and u
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
    def Jsbx(self):
        return self.__Jsbx

    @property
    def idxsbu(self):
        return self.__idxsbu

    @property
    def Jsbu(self):
        return self.__Jsbu

    # soft bounds on x at t=T
    @property
    def lsbx_e(self):
        return self.__lsbx_e

    @property
    def usbx_e(self):
        return self.__usbx_e

    @property
    def idxsbx_e(self):
        return self.__idxsbx_e

    @property
    def Jsbx_e(self):
        return self.__Jsbx_e

    # soft bounds on nonlinear constraints
    @property
    def lsh(self):
        return self.__lsh

    @property
    def ush(self):
        return self.__ush

    @property
    def idxsh(self):
        return self.__idxsh

    @property
    def Jsh(self):
        return self.__Jsh

    # soft bounds on nonlinear constraints at t=T
    @property
    def lsh_e(self):
        return self.__lsh_e

    @property
    def ush_e(self):
        return self.__ush_e

    @property
    def idxsh_e(self):
        return self.__idxsh_e


    @property
    def Jsh_e(self):
        return self.__Jsh_e

    @property
    def x0(self):
        return self.__x0

    @property
    def p(self):
        return self.__p

    # bounds on x and u
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

    @Jbx.setter
    def Jbx(self, Jbx):
        if type(Jbx) == np.ndarray:
            self.__Jbx = Jbx
        else:
            raise Exception('Invalid Jbx value. Exiting.')

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

    @Jbu.setter
    def Jbu(self, Jbu):
        if type(Jbu) == np.ndarray:
            self.__Jbu = Jbu
        else:
            raise Exception('Invalid Jbu value. Exiting.')

    # bounds on x at t=T
    @lbx_e.setter
    def lbx_e(self, lbx_e):
        if type(lbx_e) == np.ndarray:
            self.__lbx_e = lbx_e
        else:
            raise Exception('Invalid lbx_e value. Exiting.')

    @ubx_e.setter
    def ubx_e(self, ubx_e):
        if type(ubx_e) == np.ndarray:
            self.__ubx_e = ubx_e
        else:
            raise Exception('Invalid ubx_e value. Exiting.')

    @idxbx_e.setter
    def idxbx_e(self, idxbx_e):
        if type(idxbx_e) == np.ndarray:
            self.__idxbx_e = idxbx_e
        else:
            raise Exception('Invalid idxbx_e value. Exiting.')

    @Jbx_e.setter
    def Jbx_e(self, Jbx_e):
        if type(Jbx_e) == np.ndarray:
            self.__Jbx_e = Jbx_e
        else:
            raise Exception('Invalid Jbx_e value. Exiting.')

    # polytopic constraints 
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

    # polytopic constraints at t=T 
    @C_e.setter
    def C_e(self, C_e):
        if type(C_e) == np.ndarray:
            self.__C_e = C_e
        else:
            raise Exception('Invalid C_e value. Exiting.')

    @lg_e.setter
    def lg_e(self, lg_e):
        if type(lg_e) == np.ndarray:
            self.__lg_e = lg_e
        else:
            raise Exception('Invalid lg_e value. Exiting.')

    @ug_e.setter
    def ug_e(self, ug_e):
        if type(ug_e) == np.ndarray:
            self.__ug_e = ug_e
        else:
            raise Exception('Invalid ug_e value. Exiting.')

    # nonlinear constraints
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

    # nonlinear constraints at t=T
    @lh_e.setter
    def lh_e(self, lh_e):
        if type(lh_e) == np.ndarray:
            self.__lh_e = lh_e
        else:
            raise Exception('Invalid lh_e value. Exiting.')

    @uh_e.setter
    def uh_e(self, uh_e):
        if type(uh_e) == np.ndarray:
            self.__uh_e = uh_e
        else:
            raise Exception('Invalid uh_e value. Exiting.')

    # soft bounds on x and u
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

    @Jsbx.setter
    def Jsbx(self, Jsbx):
        if type(Jsbx) == np.ndarray:
            self.__Jsbx = Jsbx
        else:
            raise Exception('Invalid Jsbx value. Exiting.')


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

    @Jsbu.setter
    def Jsbu(self, Jsbu):
        if type(Jsbu) == np.ndarray:
            self.__Jsbu = Jsbu
        else:
            raise Exception('Invalid Jsbu value. Exiting.')

    # soft bounds on x at t=T
    @lsbx_e.setter
    def lsbx_e(self, lsbx_e):
        if type(lsbx_e) == np.ndarray:
            self.__lsbx_e = lsbx_e
        else:
            raise Exception('Invalid lsbx_e value. Exiting.')

    @usbx_e.setter
    def usbx_e(self, usbx_e):
        if type(usbx_e) == np.ndarray:
            self.__usbx_e = usbx_e
        else:
            raise Exception('Invalid usbx_e value. Exiting.')

    @idxsbx_e.setter
    def idxsbx_e(self, idxsbx_e):
        if type(idxsbx_e) == np.ndarray:
            self.__idxsbx_e = idxsbx_e
        else:
            raise Exception('Invalid idxsbx_e value. Exiting.')

    @Jsbx_e.setter
    def Jsbx_e(self, Jsbx_e):
        if type(Jsbx_e) == np.ndarray:
            self.__Jsbx_e = Jsbx_e
        else:
            raise Exception('Invalid Jsbx_e value. Exiting.')

    # soft bounds on nonlinear constraints
    @lsh.setter
    def lsh(self, lsh):
        if type(lsh) == np.ndarray:
            self.__lsh = lsh
        else:
            raise Exception('Invalid lsh value. Exiting.')

    @ush.setter
    def ush(self, ush):
        if type(ush) == np.ndarray:
            self.__ush = ush
        else:
            raise Exception('Invalid ush value. Exiting.')

    @idxsh.setter
    def idxsh(self, idxsh):
        if type(idxsh) == np.ndarray:
            self.__idxsh = idxsh
        else:
            raise Exception('Invalid idxsh value. Exiting.')

    # soft bounds on nonlinear constraints at t=T
    @lsh_e.setter
    def lsh_e(self, lsh_e):
        if type(lsh_e) == np.ndarray:
            self.__lsh_e = lsh_e
        else:
            raise Exception('Invalid lsh_e value. Exiting.')

    @ush_e.setter
    def ush_e(self, ush_e):
        if type(ush_e) == np.ndarray:
            self.__ush_e = ush_e
        else:
            raise Exception('Invalid ush_e value. Exiting.')

    @idxsh_e.setter
    def idxsh_e(self, idxsh_e):
        if type(idxsh_e) == np.ndarray:
            self.__idxsh_e = idxsh_e
        else:
            raise Exception('Invalid idxsh_e value. Exiting.')

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

    def set(self, attr, value):
        setattr(self, attr, value)

class ocp_nlp_solver_config:
    """
    class containing the description of the solver configuration
    """
    def __init__(self):
        self.__qp_solver        = 'PARTIAL_CONDENSING_HPIPM'  #: qp solver to be used in the NLP solver
        self.__hessian_approx   = 'GAUSS_NEWTON'              #: hessian approximation
        self.__integrator_type  = 'ERK'                       #: integrator type
        self.__tf               = None                        #: prediction horizon
        self.__nlp_solver_type  = 'SQP_RTI'                   #: NLP solver 

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

    def set(self, attr, value):
        setattr(self, attr, value)

class acados_ocp_nlp:
    """
    class containing the full description of the optimal control problem
    """
    def __init__(self):
        self.dims = ocp_nlp_dims()
        self.cost = ocp_nlp_cost()
        self.constraints = ocp_nlp_constraints()
        self.solver_config = ocp_nlp_solver_config()
        # self.con_p_name  = None 
        # self.con_p_e_name = None 
        # self.con_h_name  = None 
        # self.con_h_e_name = None 
        self.con_p   = acados_constraint() 
        self.con_p_e = acados_constraint() 
        self.con_h   = acados_constraint() 
        self.con_h_e = acados_constraint() 
        # self.constants = {}
        self.acados_include_path = []
        self.acados_lib_path = []

    def set(self, attr, value):
        # tokenize string 
        tokens = attr.split('_', 1)
        if len(tokens) > 1:
            setter_to_call = getattr(getattr(self, tokens[0]), 'set')
        else:
            setter_to_call = getattr(self, 'set')

        setter_to_call(tokens[1], value)
        return 

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
    """ Convert acados ocp nlp object to JSON format by stripping the 
    property mangling and adding array dimension info.
    ALL items of type String will be converted 
    to type ndarrray!
     
    Parameters
    ----------
    acados_ocp : class
        object of type acados_ocp_nlp.
    
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
        if (v_type == 'ndarray' or v_type__ == 'list') and (ocp_nlp_layout[k][0] != 'str'):
            dims_l = []
            dims_names = []
            dim_keys = ocp_nlp_layout[k][1]
            for item in dim_keys:
                dims_l.append(ocp_nlp_dims[item])
                dims_names.append(item)
            dims = tuple(dims_l)
            if v == []:
                # v = None
                try: 
                    v = np.reshape(v, dims)
                except:  
                    raise Exception('acados -- mismatching dimensions for field {0}. Provided data has dimensions {1}, while associated dimensions {2} are {3}'.format(out_key, [], dims_names, dims))
                # v = []
            else:
                v = np.array(v)
                v_dims = v.shape
                try: 
                    v = np.reshape(v, dims)
                except:  
                    raise Exception('acados -- mismatching dimensions for field {0}. Provided data has dimensions {1}, while associated dimensions {2} are {3}'.format(out_key, v_dims, dims_names, dims))
        out[k.replace(k, out_key)] = v
    return out

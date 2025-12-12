#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#


from acados_template import AcadosModel, AcadosOcpDims, AcadosSimDims
import numpy as np
import casadi as ca
from typing import Tuple, Union, List
import inspect, warnings

from .utils import casadi_length, is_empty, print_casadi_expression


class GnsfDims():

    def __init__(self, nx1: int, nz1: int, nuhat: int, ny: int, nout: int):

        self.__nx1 = nx1
        self.__nz1 = nz1
        self.__nuhat = nuhat
        self.__ny = ny
        self.__nout = nout

    @property
    def nx1(self):
        """GNSF: Dimension nx1."""
        return self.__nx1

    @property
    def nz1(self):
        """GNSF: Dimension nz1."""
        return self.__nz1

    @property
    def nuhat(self):
        """GNSF: Dimension nuhat."""
        return self.__nuhat

    @property
    def ny(self):
        """GNSF: Dimension ny."""
        return self.__ny

    @property
    def nout(self):
        """GNSF: Dimension nout."""
        return self.__nout

class GnsfModel():
    """
    Generalized Nonlinear Static Feedback (GNSF) model structure.

    This class represents the structured dynamic system formulation given by the
    Generalized Nonlinear Static Feedback structure (GNSF). The formulation couples:
        (i) a nonlinear static-feedback (NSF) type subsystem
            E [x1dot; z1] = A x1 + B u + C φ(L_x x1 + L_xdot x1dot + L_z z1, L_u u) + c
        (ii) a linear output system (LO)
            E_LO [x2dot; z2] = A_LO x2 + B_LO u + f_LO(x1, x1dot, z1, u) + c_LO

    Differential states are partitioned as:
        x = [x1; x2],
    algebraic variables as:
        z = [z1; z2].

    φ(·) denotes the nonlinear static-feedback function, and f_LO(·) denotes the
    nonlinearity and linear-output function evaluated in the LO subsystem.
    L_x, L_xdot, L_z, and L_u are linear input maps characterizing the affine
    argument structure of φ(·). A, B, C, E, and c are constant system matrices
    and vectors for the NSF subsystem; A_LO, E_LO, B_LO, and c_LO correspond to
    the LO subsystem.

    Parameters
    ----------
    x1 : ca.SX or ca.MX
        Differential state vector x1 ∈ R^(n_x1).
    x1dot : ca.SX or ca.MX
        Time derivative of x1.
    z1 : ca.SX or ca.MX
        Algebraic variable vector z1 ∈ R^(n_z1).
    y : ca.SX or ca.MX
        Output variable y (definition in this context: ).
    uhat : ca.SX or ca.MX
        Control variable \hat{u} = Lu * u.
    phi : ca.SX or ca.MX
        Nonlinear static-feedback function φ(·).
    f_LO : ca.SX or ca.MX
        LO nonlinear function f_LO(·).
    A : np.ndarray
        Matrix A ∈ R^{(n_x1+n_z1) x n_x1} of the NSF subsystem.
    B : np.ndarray
        Matrix B ∈ R^{(n_x1+n_z1) x n_u}.
    C : np.ndarray
        Matrix C ∈ R^{(n_x1+n_z1) x n_out}.
    E : np.ndarray
        Matrix E ∈ R^{(n_x1+n_z1) x (n_x1+n_z1)}.
    L_x : np.ndarray
        Linear input matrix L_x ∈ R^{n_y x n_x1}.
    L_xdot : np.ndarray
        Linear input matrix L_xdot ∈ R^{n_y x n_x1}.
    L_u : np.ndarray
        Linear input matrix L_u ∈ R^{n_y x n_u}.
    L_z : np.ndarray
        Linear input matrix L_z ∈ R^{n_y x n_z1}.
    A_LO : np.ndarray
        Matrix A_LO of the LO subsystem.
    c : np.ndarray
        Constant vector c ∈ R^{n_x1+n_z1}.
    E_LO : np.ndarray
        Matrix E_LO of the LO subsystem.
    B_LO : np.ndarray
        Matrix B_LO of the LO subsystem.
    c_LO : np.ndarray
        Constant vector c_LO.
    idx_perm_x : np.ndarray
        Permutation/pivot indices for x-partitioning.
    idx_perm_z : np.ndarray
        Permutation/pivot indices for z-partitioning.
    purely_linear : bool
        Indicates whether the model contains no nonlinear φ(·) component.
    nontrivial_f_LO : bool
        Indicates whether f_LO contains a nontrivial nonlinear contribution.
    """

    def __init__(self,
                 x: Union[ca.SX, ca.MX],
                 u: Union[ca.SX, ca.MX],
                 z: Union[ca.SX, ca.MX],
                 xdot: Union[ca.SX, ca.MX],
                 p: Union[ca.SX, ca.MX],
                 p_global: Union[ca.SX, ca.MX],
                 y: Union[ca.SX, ca.MX],
                 uhat: Union[ca.SX, ca.MX],
                 phi: Union[ca.SX, ca.MX],
                 f_LO: Union[ca.SX, ca.MX],
                 A: np.ndarray,
                 B: np.ndarray,
                 C: np.ndarray,
                 E: np.ndarray,
                 A_LO: np.ndarray,
                 c: np.ndarray,
                 E_LO: np.ndarray,
                 B_LO: np.ndarray,
                 c_LO: np.ndarray,
                 idx_perm_x: np.ndarray,
                 idx_perm_z: np.ndarray,
                 ):

        # symbolics
        self.__x = x
        self.__u = u
        self.__z = z
        self.__xdot = xdot
        self.__p = p
        self.__p_global = p_global

        # expressions
        self.__y = y
        self.__uhat = uhat
        self.__phi = phi
        self.__f_LO = f_LO

        # matrices and vectors
        self.__A = A
        self.__B = B
        self.__C = C
        self.__E = E
        self.__A_LO = A_LO
        self.__c = c
        self.__E_LO = E_LO
        self.__B_LO = B_LO
        self.__c_LO = c_LO

        self.__idx_perm_x = idx_perm_x
        self.__idx_perm_z = idx_perm_z

        # everything below will be detected
        self.__x1 = None
        self.__z1 = None
        self.__x1dot = None

        self.__L_x = None
        self.__L_xdot = None
        self.__L_u = None
        self.__L_z = None

        self.__purely_linear = None
        self.__nontrivial_f_LO = None

        self.__dims = None

        self._detect_dims()
        self._make_consistent()

    @property
    def x(self):
        """GNSF: Symbolic state vector x."""
        return self.__x

    @property
    def u(self):
        """GNSF: Symbolic control vector u."""
        return self.__u

    @property
    def z(self):
        """GNSF: Symbolic algebraic variable vector z."""
        return self.__z

    @property
    def xdot(self):
        """GNSF: Symbolic state derivative vector xdot."""
        return self.__xdot

    @property
    def p(self):
        """GNSF: Symbolic parameter vector p."""
        return self.__p

    @property
    def p_global(self):
        """GNSF: Symbolic global parameter vector p_global."""
        return self.__p_global

    @property
    def y(self):
        """GNSF: Symbolic expression for y in GNSF formulation."""
        return self.__y

    @property
    def uhat(self):
        """GNSF: Symbolic expression for uhat in GNSF formulation."""
        return self.__uhat

    @property
    def x1(self):
        """GNSF: Symbolic expression for x1 in GNSF formulation."""
        return self.__x1

    @property
    def z1(self):
        """GNSF: Symbolic expression for z1 in GNSF formulation."""
        return self.__z1

    @property
    def x1dot(self):
        """GNSF: Symbolic expression for x1dot in GNSF formulation."""
        return self.__x1dot

    @property
    def phi(self):
        """GNSF: Symbolic expression for phi in GNSF formulation."""
        return self.__phi

    @property
    def f_LO(self):
        """GNSF: Symbolic expression for f_LO (linear output function) in GNSF formulation."""
        return self.__f_LO

    @property
    def A(self):
        """GNSF: Matrix A in GNSF formulation."""
        return self.__A

    @property
    def B(self):
        """GNSF: Matrix B in GNSF formulation."""
        return self.__B

    @property
    def C(self):
        """GNSF: Matrix C in GNSF formulation."""
        return self.__C

    @property
    def E(self):
        """GNSF: Matrix E in GNSF formulation."""
        return self.__E

    @property
    def L_x(self):
        """GNSF: Matrix L_x in GNSF formulation."""
        return self.__L_x

    @property
    def L_xdot(self):
        """GNSF: Matrix L_xdot in GNSF formulation."""
        return self.__L_xdot

    @property
    def L_u(self):
        """GNSF: Matrix L_u in GNSF formulation."""
        return self.__L_u

    @property
    def L_z(self):
        """GNSF: Matrix L_z in GNSF formulation."""
        return self.__L_z

    @property
    def A_LO(self):
        """GNSF: Matrix A_LO (linear output) in GNSF formulation."""
        return self.__A_LO

    @property
    def c(self):
        """GNSF: Vector c in GNSF formulation."""
        return self.__c

    @property
    def E_LO(self):
        """GNSF: Matrix E_LO (linear output) in GNSF formulation."""
        return self.__E_LO

    @property
    def B_LO(self):
        """GNSF: Matrix B_LO (linear output) in GNSF formulation."""
        return self.__B_LO

    @property
    def nontrivial_f_LO(self):
        """GNSF: Flag indicating whether GNSF structure has nontrivial f_LO."""
        return self.__nontrivial_f_LO

    @property
    def purely_linear(self):
        """GNSF: Flag indicating whether GNSF structure is purely linear."""
        return self.__purely_linear

    @property
    def idx_perm_x(self):
        """GNSF: Permutation index vector for x in GNSF formulation."""
        return self.__idx_perm_x

    @property
    def idx_perm_z(self):
        """GNSF: Permutation index vector for x in GNSF formulation."""
        return self.__idx_perm_z

    @property
    def c_LO(self):
        """GNSF: Vector c_LO (linear output) in GNSF formulation."""
        return self.__c_LO

    @property
    def dims(self):
        """GNSF: Dimensions of GNSF model."""
        return self.__dims


    def _detect_dims(self,):
        """
        Detect dimensions of GNSF model.
        """

        nx1 = self.A.shape[1]
        nz1 = self.A.shape[0] - nx1
        nuhat = self.uhat.rows()
        ny = self.y.rows()
        nout = self.C.shape[1]

        self.__dims = GnsfDims(nx1, nz1, nuhat, ny, nout)


    def _make_consistent(self,):
        """
        Make sure the GNSF model is consistent.
        """

        # detect linear input system
        if not (ca.is_linear(self.y, self.x) and ca.is_linear(self.y, self.xdot) and ca.is_linear(self.y, self.z)):
            raise ValueError("y must be linear in x, xdot, z")

        if not (ca.is_linear(self.uhat, self.u)):
            raise ValueError("uhat must be linear in u")

        self.__L_x = ca.evalf(ca.jacobian(self.y, self.x)).full()
        self.__L_xdot = ca.evalf(ca.jacobian(self.y, self.xdot)).full()

        if not is_empty(self.z):
            self.__L_z = ca.evalf(ca.jacobian(self.y, self.z)).full()
        else:
            self.__L_z = np.zeros((self.dims.ny, 0))
        self.__L_u = ca.evalf(ca.jacobian(self.uhat, self.u)).full()


        # detect flags
        self.__nontrivial_f_LO = not is_empty(self.f_LO) and not self.f_LO.is_zero()
        self.__purely_linear = self.dims.nx1 == 0 and self.dims.nz1 == 0 and not self.nontrivial_f_LO

        breakpoint()
        # TODO: sanity checks
        self.__ipipv_x = idx_perm_to_ipiv(self.idx_perm_x)
        self.__ipipv_z = idx_perm_to_ipiv(self.idx_perm_z)

        pass


    @classmethod
    def detect(cls, model: AcadosModel, dims: Union[AcadosOcpDims, AcadosSimDims]):
        """
        Detect GnsfModel from AcadosModel.
        """

        gnsf_model = cls()

        return gnsf_model


    def to_dict(self):
        """
        Convert GnsfModel to dictionary.

        Returns:
            dict: dictionary representation of GnsfModel
        """

        model_dict = {}

        for k, _ in inspect.getmembers(type(self), lambda v: isinstance(v, property)):
            v = getattr(self, k)
            if isinstance(v, (ca.SX, ca.MX)):
                model_dict[k] = repr(v) # only for debugging
            elif isinstance(v, GnsfDims):
                model_dict[k] = dict(getattr(self, k).__dict__)
            else:
                model_dict[k] = v

        model_dict['serialized_expressions'], model_dict['expression_names'] = self.serialize()

        return model_dict


    def serialize(self) -> Tuple[str, List[str]]:
        """
        Serialize the CasADi expressions.
        """

        serializer = ca.StringSerializer()
        expression_names = []

        for k, _ in inspect.getmembers(type(self), lambda v: isinstance(v, property)):
            v = getattr(self, k)
            if isinstance(v, (ca.SX, ca.MX)):
                serializer.pack(v)
                expression_names.append(k)

        return serializer.encode(), expression_names


    def deserialize(self, s: str, expression_names: list) -> None:
        """
        Deserialize the CasADi expressions.
        """

        deserializer = ca.StringDeserializer(s)

        for name in expression_names:
            setattr(self, name, deserializer.unpack())


    @classmethod
    def from_dict(cls, model_dict: dict) -> 'GnsfModel':
        """
        Create an GnsfModel from a dictionary.
        Values that correspond to the empty list are ignored.
        """

        model = cls()

        expression_names =  model_dict.get('expression_names')
        serialized_expressions = model_dict.get('serialized_expressions')

        if expression_names is None or serialized_expressions is None:
            raise ValueError("Dictionary does not contain serialized expressions.")

        # loop over all properties
        for attr, _ in inspect.getmembers(type(model), lambda v: isinstance(v, property)):

            value = model_dict.get(attr)

            # expressions are expected to be None
            if value is None and attr not in expression_names:
                warnings.warn(f"Attribute {attr} not in dictionary.")
            else:
                try:
                    # check whether value is not the empty list and not a CasADi symbol/expression
                    if not (isinstance(value, list) and not value) and not attr in expression_names:
                        setattr(model, attr, value)
                except Exception as e:
                    Exception("Failed to load attribute {attr} from dictionary:\n" + repr(e))

        model.deserialize(model_dict['serialized_expressions'], model_dict['expression_names'])

        return model



def check_reformulation(model, gnsf, print_info):

    ## Description:
    # this function takes the implicit ODE/ index-1 DAE and a gnsf structure
    # to evaluate both models at num_eval random points x0, x0dot, z0, u0
    # if for all points the relative error is <= TOL, the function will return::
    # 1, otherwise it will give an error.

    TOL = 1e-14
    num_eval = 10

    # get dimensions
    nx = gnsf["nx"]
    nu = gnsf["nu"]
    nz = gnsf["nz"]
    nx1 = gnsf["nx1"]
    nx2 = gnsf["nx2"]
    nz1 = gnsf["nz1"]
    nz2 = gnsf["nz2"]
    n_out = gnsf["n_out"]

    # get model matrices
    A = gnsf["A"]
    B = gnsf["B"]
    C = gnsf["C"]
    E = gnsf["E"]
    c = gnsf["c"]

    L_x = gnsf["L_x"]
    L_xdot = gnsf["L_xdot"]
    L_z = gnsf["L_z"]
    L_u = gnsf["L_u"]

    A_LO = gnsf["A_LO"]
    E_LO = gnsf["E_LO"]
    B_LO = gnsf["B_LO"]
    c_LO = gnsf["c_LO"]

    I_x1 = range(nx1)
    I_x2 = range(nx1, nx)

    I_z1 = range(nz1)
    I_z2 = range(nz1, nz)

    idx_perm_f = gnsf["idx_perm_f"]

    # get casadi variables
    x = gnsf["x"]
    xdot = gnsf["xdot"]
    z = gnsf["z"]
    u = gnsf["u"]
    y = gnsf["y"]
    uhat = gnsf["uhat"]
    p = gnsf["p"]

    # create functions
    impl_dae_fun = ca.Function("impl_dae_fun", [x, xdot, u, z, p], [model.f_impl_expr])
    phi_fun = ca.Function("phi_fun", [y, uhat, p], [gnsf["phi_expr"]])
    f_lo_fun = ca.Function(
        "f_lo_fun", [x[range(nx1)], xdot[range(nx1)], z, u, p], [gnsf["f_lo_expr"]]
    )

    for i_check in range(num_eval):
        # generate random values
        x0 = np.random.rand(nx, 1)
        x0dot = np.random.rand(nx, 1)
        z0 = np.random.rand(nz, 1)
        u0 = np.random.rand(nu, 1)

        if gnsf["ny"] > 0:
            y0 = L_x @ x0[I_x1] + L_xdot @ x0dot[I_x1] + L_z @ z0[I_z1]
        else:
            y0 = []
        if gnsf["nuhat"] > 0:
            uhat0 = L_u @ u0
        else:
            uhat0 = []

        # eval functions
        p0 = np.random.rand(gnsf["np"], 1)
        f_impl_val = impl_dae_fun(x0, x0dot, u0, z0, p0).full()
        phi_val = phi_fun(y0, uhat0, p0)
        f_lo_val = f_lo_fun(x0[I_x1], x0dot[I_x1], z0[I_z1], u0, p0)

        f_impl_val = f_impl_val[idx_perm_f]
        # eval gnsf
        if n_out > 0:
            C_phi = C @ phi_val
        else:
            C_phi = np.zeros((nx1 + nz1, 1))
        try:
            gnsf_val1 = (
                A @ x0[I_x1] + B @ u0 + C_phi + c - E @ ca.vertcat(x0dot[I_x1], z0[I_z1])
            )
            # gnsf_1 = (A @ x[I_x1] + B @ u + C_phi + c - E @ ca.vertcat(xdot[I_x1], z[I_z1]))
        except:
            breakpoint()

        if nx2 > 0:  # eval LOS:
            gnsf_val2 = (
                A_LO @ x0[I_x2]
                + B_LO @ u0
                + c_LO
                + f_lo_val
                - E_LO @ ca.vertcat(x0dot[I_x2], z0[I_z2])
            )
            gnsf_val = ca.vertcat(gnsf_val1, gnsf_val2).full()
        else:
            gnsf_val = gnsf_val1.full()
        # compute error and check
        rel_error = np.linalg.norm(f_impl_val - gnsf_val) / np.linalg.norm(f_impl_val)

        if rel_error > TOL:
            print("transcription failed rel_error > TOL")
            abs_error = gnsf_val - f_impl_val
            # T = table(f_impl_val, gnsf_val, abs_error)
            # print(T)
            print("abs_error:", abs_error)
            #         error('transcription failed rel_error > TOL')
            #         check = 0
            breakpoint()

    if print_info:
        print(" ")
        print("model reformulation checked: relative error <= TOL = ", str(TOL))
        print(" ")
        check = 1
    ## helpful for debugging:
    # # use in calling function and compare
    # # compare f_impl(i) with gnsf_val1(i)
    #

    #     nx  = gnsf['nx']
    #     nu  = gnsf['nu']
    #     nz  = gnsf['nz']
    #     nx1 = gnsf['nx1']
    #     nx2 = gnsf['nx2']
    #
    #         A  = gnsf['A']
    #     B  = gnsf['B']
    #     C  = gnsf['C']
    #     E  = gnsf['E']
    #     c  = gnsf['c']
    #
    #     L_x    = gnsf['L_x']
    #     L_z    = gnsf['L_z']
    #     L_xdot = gnsf['L_xdot']
    #     L_u    = gnsf['L_u']
    #
    #     A_LO = gnsf['A_LO']
    #
    #     x0 = rand(nx, 1)
    #     x0dot = rand(nx, 1)
    #     z0 = rand(nz, 1)
    #     u0 = rand(nu, 1)
    #     I_x1 = range(nx1)
    #     I_x2 = nx1+range(nx)
    #
    #     y0 = L_x @ x0[I_x1] + L_xdot @ x0dot[I_x1] + L_z @ z0
    #     uhat0 = L_u @ u0
    #
    #     gnsf_val1 = (A @ x[I_x1] + B @ u + #         C @ phi_current + c) - E @ [xdot[I_x1] z]
    #     gnsf_val1 = gnsf_val1.ca.simplify()
    #
    # #     gnsf_val2 = A_LO @ x[I_x2] + gnsf['f_lo_fun'](x[I_x1], xdot[I_x1], z, u) - xdot[I_x2]
    #     gnsf_val2 =  A_LO @ x[I_x2] + gnsf['f_lo_fun'](x[I_x1], xdot[I_x1], z, u) - xdot[I_x2]
    #
    #
    #     gnsf_val = [gnsf_val1 gnsf_val2]
    #     gnsf_val = gnsf_val.ca.simplify()
    #     dyn_expr_f = dyn_expr_f.ca.simplify()

    return check




def detect_gnsf_structure(model: AcadosModel, dims: Union[AcadosSimDims, AcadosOcpDims],
                        print_info: bool = True,
                        detect_LOS: bool = True,
                        check_E_invertibility: bool = True,
                    ):

    ## Description
    # This function takes a CasADi implicit ODE or index-1 DAE model "model"
    # consisting of a CasADi expression f_impl in the symbolic CasADi
    # variables x, xdot, u, z, (and possibly parameters p), which are also part
    # of the model, as well as a model name.
    # It will create a struct "gnsf" containing all information needed to use
    # it with the gnsf integrator in acados.
    # Additionally it will create the struct "reordered_model" which contains
    # the permuted state vector and permuted f_impl, in which additionally some
    # functions, which were made part of the linear output system of the gnsf,
    # have changed signs.

    #   print_info: if extensive information on how the model is processed
    #       is printed to the console.
    #   check_E_invertibility: if the transcription method should check if the
    #       assumption that the main blocks of the matrix gnsf.E are invertible
    #       holds. If not, the method will try to reformulate the gnsf model
    #       with a different model, such that the assumption holds.

    # acados_root_dir = getenv('ACADOS_INSTALL_DIR')

    if not is_empty(model.p_global) and ca.depends_on(model.f_impl_expr, model.p_global):
        NotImplementedError("GNSF does not support global parameters")

    ## Reformulate implicit index-1 DAE into GNSF form
    # (Generalized nonlinear static feedback)
    gnsf = determine_trivial_gnsf_transcription(model, dims, print_info)
    gnsf = detect_affine_terms_reduce_nonlinearity(gnsf, model, print_info)

    if detect_LOS:
        gnsf = reformulate_with_LOS(model, gnsf, print_info)

    if check_E_invertibility:
        gnsf = reformulate_with_invertible_E_mat(gnsf, model, print_info)

    # detect purely linear model
    if gnsf["nx1"] == 0 and gnsf["nz1"] == 0 and gnsf["nontrivial_f_LO"] == 0:
        gnsf["purely_linear"] = 1
    else:
        gnsf["purely_linear"] = 0

    structure_detection_print_summary(gnsf, model)
    check_reformulation(model, gnsf, print_info)

    ## copy relevant fields from gnsf to model
    model_name = model.name

    phi = gnsf["phi_expr"]
    y = gnsf["y"]
    uhat = gnsf["uhat"]
    p = gnsf["p"]

    jac_phi_y = ca.jacobian(phi, y)
    jac_phi_uhat = ca.jacobian(phi, uhat)

    phi_fun = ca.Function(f"{model_name}_gnsf_phi_fun", [y, uhat, p], [phi])
    model.phi_fun = phi_fun
    model.phi_fun_jac_y = ca.Function(
        f"{model_name}_gnsf_phi_fun_jac_y", [y, uhat, p], [phi, jac_phi_y]
    )
    model.phi_jac_y_uhat = ca.Function(
        f"{model_name}_gnsf_phi_jac_y_uhat", [y, uhat, p], [jac_phi_y, jac_phi_uhat]
    )

    x1 = model.x[gnsf["idx_perm_x"][: gnsf["nx1"]]]
    x1dot = model.xdot[gnsf["idx_perm_x"][: gnsf["nx1"]]]
    if gnsf["nz1"] > 0:
        z1 = model.z[gnsf["idx_perm_z"][: gnsf["nz1"]]]
    else:
        z1 = ca.SX.sym("z1", 0, 0)

    # flags
    model.gnsf_nontrivial_f_LO = gnsf['nontrivial_f_LO']
    model.gnsf_purely_linear = gnsf['purely_linear']

    model.gnsf_model = GnsfModel(
        y=gnsf["y"],
        uhat=gnsf["uhat"],
        x1=x1,
        z1=z1,
        x1dot=x1dot,
        phi=gnsf["phi_expr"],
        f_LO=gnsf["f_lo_expr"],
        A=gnsf["A"],
        B=gnsf["B"],
        C=gnsf["C"],
        E=gnsf["E"],
        L_x=gnsf["L_x"],
        L_xdot=gnsf["L_xdot"],
        L_u=gnsf["L_u"],
        L_z=gnsf["L_z"],
        A_LO=gnsf["A_LO"],
        c=gnsf["c"],
        E_LO=gnsf["E_LO"],
        B_LO=gnsf["B_LO"],
        c_LO=gnsf["c_LO"],
        idx_perm_x=gnsf["idx_perm_x"],
        idx_perm_z=gnsf["idx_perm_z"],
        purely_linear=gnsf["purely_linear"],
        nontrivial_f_LO=gnsf["nontrivial_f_LO"],
    )

    return



def detect_affine_terms_reduce_nonlinearity(gnsf, model: AcadosModel, print_info):

    ## Description
    # this function takes a gnsf structure with trivial model matrices (A, B,
    # E, c are zeros, and C is eye).
    # It detects all affine linear terms and sets up an equivalent model in the
    # GNSF structure, where all affine linear terms are modeled through the
    # matrices A, B, E, c and the linear output system (LOS) is empty.
    # NOTE: model is just taken as an argument to check equivalence of the
    # models within the function.

    if print_info:
        print(" ")
        print("====================================================================")
        print(" ")
        print("============  Detect affine-linear dependencies   ==================")
        print(" ")
        print("====================================================================")
        print(" ")
    # symbolics
    x = gnsf["x"]
    xdot = gnsf["xdot"]
    u = gnsf["u"]
    z = gnsf["z"]

    # dimensions
    nx = gnsf["nx"]
    nu = gnsf["nu"]
    nz = gnsf["nz"]

    ny_old = gnsf["ny"]
    nuhat_old = gnsf["nuhat"]

    ## Represent all affine dependencies through the model matrices A, B, E, c
    ## determine A
    n_nodes_current = ca.n_nodes(gnsf["phi_expr"])

    for ii in range(casadi_length(gnsf["phi_expr"])):
        fii = gnsf["phi_expr"][ii]
        for ix in range(nx):
            var = x[ix]
            varname = var.name()
            # symbolic jacobian of fii w.r.t. xi
            jac_fii_xi = ca.jacobian(fii, var)
            if jac_fii_xi.is_constant():
                # jacobian value
                jac_fii_xi_fun = ca.Function("jac_fii_xi_fun", [x[1]], [jac_fii_xi])
                # x[1] as input just to have a scalar input and call the function as follows:
                gnsf["A"][ii, ix] = jac_fii_xi_fun(0).full()
            else:
                gnsf["A"][ii, ix] = 0
                if print_info:
                    print(f"phi({ii}) is nonlinear in x({ix}) = {varname}")
                    print(fii)
                    print("-----------------------------------------------------")
    f_next = gnsf["phi_expr"] - gnsf["A"] @ x
    f_next = ca.simplify(f_next)
    n_nodes_next = ca.n_nodes(f_next)

    if print_info:
        print("\n")
        print(f"determined matrix A:")
        print(gnsf["A"])
        print(f"reduced nonlinearity from  {n_nodes_current} to {n_nodes_next} nodes")
    # assert(n_nodes_current >= n_nodes_next,'n_nodes_current >= n_nodes_next FAILED')
    gnsf["phi_expr"] = f_next

    check_reformulation(model, gnsf, print_info)

    ## determine B
    n_nodes_current = ca.n_nodes(gnsf["phi_expr"])

    for ii in range(casadi_length(gnsf["phi_expr"])):
        fii = gnsf["phi_expr"][ii]
        for iu in range(nu):
            var = u[iu]
            varname = var.name()
            # symbolic jacobian of fii w.r.t. ui
            jac_fii_ui = ca.jacobian(fii, var)
            if jac_fii_ui.is_constant():  # i.e. hessian is structural zero:
                # jacobian value
                jac_fii_ui_fun = ca.Function("jac_fii_ui_fun", [x[1]], [jac_fii_ui])
                gnsf["B"][ii, iu] = jac_fii_ui_fun(0).full()
            else:
                gnsf["B"][ii, iu] = 0
                if print_info:
                    print(f"phi({ii}) is nonlinear in u(", str(iu), ") = ", varname)
                    print(fii)
                    print("-----------------------------------------------------")
    f_next = gnsf["phi_expr"] - gnsf["B"] @ u
    f_next = ca.simplify(f_next)
    n_nodes_next = ca.n_nodes(f_next)

    if print_info:
        print("\n")
        print(f"determined matrix B:")
        print(gnsf["B"])
        print(f"reduced nonlinearity from  {n_nodes_current} to {n_nodes_next} nodes")

    gnsf["phi_expr"] = f_next

    check_reformulation(model, gnsf, print_info)

    ## determine E
    n_nodes_current = ca.n_nodes(gnsf["phi_expr"])
    k = ca.vertcat(xdot, z)

    for ii in range(casadi_length(gnsf["phi_expr"])):
        fii = gnsf["phi_expr"][ii]
        for ik in range(casadi_length(k)):
            # symbolic jacobian of fii w.r.t. ui
            var = k[ik]
            varname = var.name
            jac_fii_ki = ca.jacobian(fii, var)
            if jac_fii_ki.is_constant():
                # jacobian value
                jac_fii_ki_fun = ca.Function("jac_fii_ki_fun", [x[1]], [jac_fii_ki])
                gnsf["E"][ii, ik] = -jac_fii_ki_fun(0).full()
            else:
                gnsf["E"][ii, ik] = 0
                if print_info:
                    print(f"phi( {ii}) is nonlinear in xdot_z({ik}) = ", varname)
                    print(fii)
                    print("-----------------------------------------------------")
    f_next = gnsf["phi_expr"] + gnsf["E"] @ k
    f_next = ca.simplify(f_next)
    n_nodes_next = ca.n_nodes(f_next)

    if print_info:
        print("\n")
        print(f"determined matrix E:")
        print(gnsf["E"])
        print(f"reduced nonlinearity from {n_nodes_current} to {n_nodes_next} nodes")

    gnsf["phi_expr"] = f_next
    check_reformulation(model, gnsf, print_info)

    ## determine constant term c

    n_nodes_current = ca.n_nodes(gnsf["phi_expr"])
    for ii in range(casadi_length(gnsf["phi_expr"])):
        fii = gnsf["phi_expr"][ii]
        if fii.is_constant():
            # function value goes into c
            fii_fun = ca.Function("fii_fun", [x[1]], [fii])
            gnsf["c"][ii] = fii_fun(0).full()
        else:
            gnsf["c"][ii] = 0
            if print_info:
                print(f"phi(", str(ii), ") is NOT constant")
                print(fii)
                print("-----------------------------------------------------")
    gnsf["phi_expr"] = gnsf["phi_expr"] - gnsf["c"]
    gnsf["phi_expr"] = ca.simplify(gnsf["phi_expr"])
    n_nodes_next = ca.n_nodes(gnsf["phi_expr"])

    if print_info:
        print("\n")
        print(f"determined vector c:")
        print(gnsf["c"])
        print(f"reduced nonlinearity from {n_nodes_current} to {n_nodes_next} nodes")

    check_reformulation(model, gnsf, print_info)

    ## determine nonlinearity & corresponding matrix C
    ## Reduce dimension of phi
    n_nodes_current = ca.n_nodes(gnsf["phi_expr"])
    ind_non_zero = []
    for ii in range(casadi_length(gnsf["phi_expr"])):
        fii = gnsf["phi_expr"][ii]
        fii = ca.simplify(fii)
        if not fii.is_zero():
            ind_non_zero = list(set.union(set(ind_non_zero), set([ii])))
    gnsf["phi_expr"] = gnsf["phi_expr"][ind_non_zero]

    # C
    gnsf["C"] = np.zeros((nx + nz, len(ind_non_zero)))
    for ii in range(len(ind_non_zero)):
        gnsf["C"][ind_non_zero[ii], ii] = 1
    gnsf = determine_input_nonlinearity_function(gnsf)
    n_nodes_next = ca.n_nodes(gnsf["phi_expr"])

    if print_info:
        print(" ")
        print("determined matrix C:")
        print(gnsf["C"])
        print(
            "---------------------------------------------------------------------------------"
        )
        print(
            "------------- Success: Affine linear terms detected -----------------------------"
        )
        print(
            "---------------------------------------------------------------------------------"
        )
        print(
            f'reduced nonlinearity dimension n_out from  {nx+nz}  to  {gnsf["n_out"]}'
        )
        print(f"reduced nonlinearity from  {n_nodes_current} to {n_nodes_next} nodes")
        print(" ")
        print("phi now reads as:")
        print_casadi_expression(gnsf["phi_expr"])

    ## determine input of nonlinearity function
    check_reformulation(model, gnsf, print_info)

    gnsf["ny"] = casadi_length(gnsf["y"])
    gnsf["nuhat"] = casadi_length(gnsf["uhat"])

    if print_info:
        print(
            "-----------------------------------------------------------------------------------"
        )
        print(" ")
        print(
            f"reduced input ny    of phi from  ",
            str(ny_old),
            "   to  ",
            str(gnsf["ny"]),
        )
        print(
            f"reduced input nuhat of phi from  ",
            str(nuhat_old),
            "   to  ",
            str(gnsf["nuhat"]),
        )
        print(
            "-----------------------------------------------------------------------------------"
        )

    # if print_info:
    #     print(f"gnsf: {gnsf}")

    return gnsf




def determine_input_nonlinearity_function(gnsf):

    ## Description
    # this function takes a structure gnsf and updates the matrices L_x,
    # L_xdot, L_z, L_u and CasADi vectors y, uhat of this structure as follows:

    # given a CasADi expression phi_expr, which may depend on the variables
    # (x1, x1dot, z, u), this function determines a vector y (uhat) consisting
    # of all components of (x1, x1dot, z) (respectively u) that enter phi_expr.
    # Additionally matrices L_x, L_xdot, L_z, L_u are determined such that
    #           y    = L_x * x + L_xdot * xdot + L_z * z
    #           uhat = L_u * u
    # Furthermore the dimensions ny, nuhat, n_out are updated

    ## y
    y = ca.SX.sym('y', 0, 0)
    # components of x1
    for ii in range(gnsf["nx1"]):
        if ca.which_depends(gnsf["phi_expr"], gnsf["x"][ii])[0]:
            y = ca.vertcat(y, gnsf["x"][ii])
        # else:
        # x[ii] is not part of y
    # components of x1dot
    for ii in range(gnsf["nx1"]):
        if ca.which_depends(gnsf["phi_expr"], gnsf["xdot"][ii])[0]:
            print(gnsf["phi_expr"], "depends on", gnsf["xdot"][ii])
            y = ca.vertcat(y, gnsf["xdot"][ii])
        # else:
        # xdot[ii] is not part of y
    # components of z
    for ii in range(gnsf["nz1"]):
        if ca.which_depends(gnsf["phi_expr"], gnsf["z"][ii])[0]:
            y = ca.vertcat(y, gnsf["z"][ii])
        # else:
        # z[ii] is not part of y
    ## uhat
    uhat = ca.SX.sym('uhat', 0, 0)
    # components of u
    for ii in range(gnsf["nu"]):
        if ca.which_depends(gnsf["phi_expr"], gnsf["u"][ii])[0]:
            uhat = ca.vertcat(uhat, gnsf["u"][ii])
        # else:
        # u[ii] is not part of uhat
    ## generate gnsf['phi_expr_fun']
    # linear input matrices
    if is_empty(y):
        gnsf["L_x"] = []
        gnsf["L_xdot"] = []
        gnsf["L_u"] = []
        gnsf["L_z"] = []
    else:
        dummy = ca.SX.sym("dummy_input", 0)
        L_x_fun = ca.Function(
            "L_x_fun", [dummy], [ca.jacobian(y, gnsf["x"][range(gnsf["nx1"])])]
        )
        L_xdot_fun = ca.Function(
            "L_xdot_fun", [dummy], [ca.jacobian(y, gnsf["xdot"][range(gnsf["nx1"])])]
        )
        L_z_fun = ca.Function(
            "L_z_fun", [dummy], [ca.jacobian(y, gnsf["z"][range(gnsf["nz1"])])]
        )
        L_u_fun = ca.Function("L_u_fun", [dummy], [ca.jacobian(uhat, gnsf["u"])])

        gnsf["L_x"] = L_x_fun(0).full()
        gnsf["L_xdot"] = L_xdot_fun(0).full()
        gnsf["L_u"] = L_u_fun(0).full()
        gnsf["L_z"] = L_z_fun(0).full()
    gnsf["y"] = y
    gnsf["uhat"] = uhat

    gnsf["ny"] = casadi_length(y)
    gnsf["nuhat"] = casadi_length(uhat)
    gnsf["n_out"] = casadi_length(gnsf["phi_expr"])

    return gnsf




def determine_trivial_gnsf_transcription(model: AcadosModel, dims: Union[AcadosSimDims, AcadosOcpDims], print_info):
    ## Description
    # this function takes a model of an implicit ODE/ index-1 DAE and sets up
    # an equivalent model in the GNSF structure, with empty linear output
    # system and trivial model matrices, i.e. A, B, E, c are zeros, and C is
    # eye. - no structure is exploited

    # initial print
    print("*****************************************************************")
    print(" ")
    print(f"******      Restructuring {model.name} model    ***********")
    print(" ")
    print("*****************************************************************")

    # load model
    f_impl_expr = model.f_impl_expr

    model_name_prefix = model.name

    # x
    x = model.x
    nx = dims.nx
    # check type
    if isinstance(x[0], ca.SX):
        isSX = True
    else:
        print("GNSF detection only works for SX CasADi type!!!")
        breakpoint()
    # xdot
    xdot = model.xdot
    # u
    nu = dims.nu
    if nu == 0:
        u = ca.SX.sym("u", 0, 0)
    else:
        u = model.u

    nz = dims.nz
    if nz == 0:
        z = ca.SX.sym("z", 0, 0)
    else:
        z = model.z

    p = model.p
    nparam = dims.np

    # avoid SX of size 0x1
    if casadi_length(u) == 0:
        u = ca.SX.sym("u", 0, 0)
        nu = 0
    ## initialize gnsf struct
    # dimensions
    gnsf = {"nx": nx, "nu": nu, "nz": nz, "np": nparam}
    gnsf["nx1"] = nx
    gnsf["nx2"] = 0
    gnsf["nz1"] = nz
    gnsf["nz2"] = 0
    gnsf["nuhat"] = nu
    gnsf["ny"] = 2 * nx + nz

    gnsf["phi_expr"] = f_impl_expr
    gnsf["A"] = np.zeros((nx + nz, nx))
    gnsf["B"] = np.zeros((nx + nz, nu))
    gnsf["E"] = np.zeros((nx + nz, nx + nz))
    gnsf["c"] = np.zeros((nx + nz, 1))
    gnsf["C"] = np.eye(nx + nz)
    gnsf["name"] = model_name_prefix

    gnsf["x"] = x
    gnsf["xdot"] = xdot
    gnsf["z"] = z
    gnsf["u"] = u
    gnsf["p"] = p

    gnsf = determine_input_nonlinearity_function(gnsf)

    gnsf["A_LO"] = []
    gnsf["E_LO"] = []
    gnsf["B_LO"] = []
    gnsf["c_LO"] = []
    gnsf["f_lo_expr"] = []

    # permutation
    gnsf["idx_perm_x"] = range(nx)  # matlab-style)
    gnsf["idx_perm_z"] = range(nz)
    gnsf["idx_perm_f"] = range((nx + nz))

    gnsf["nontrivial_f_LO"] = 0

    check_reformulation(model, gnsf, print_info)
    if print_info:
        print(f"Success: Set up equivalent GNSF model with trivial matrices")
        print(" ")
    if print_info:
        print("-----------------------------------------------------------")
        print(" ")
        print(f"reduced input ny    of phi from  {2 * nx + nz}   to  {gnsf['ny']}")
        print(f"reduced input nuhat of phi from  {nu}   to  {gnsf['nuhat']}")
        print("-----------------------------------------------------------")
    return gnsf




def reformulate_with_invertible_E_mat(gnsf, model, print_info):
    ## Description
    # this function checks that the necessary condition to apply the gnsf
    # structure exploiting integrator to a model, namely that the matrices E11,
    # E22 are invertible holds.
    # if this is not the case, it will make these matrices invertible and add:
    # corresponding terms, to the term C * phi, such that the obtained model is
    # still equivalent

    # check invertibility of E11, E22 and reformulate if needed:
    ind_11 = range(gnsf["nx1"])
    ind_22 = range(gnsf["nx1"], gnsf["nx1"] + gnsf["nz1"])

    if print_info:
        print(" ")
        print("----------------------------------------------------")
        print("checking rank of E11 and E22")
        print("----------------------------------------------------")
    ## check if E11, E22 are invertible:
    z_check = False
    if gnsf["nz1"] > 0:
        z_check = (
            np.linalg.matrix_rank(gnsf["E"][np.ix_(ind_22, ind_22)]) != gnsf["nz1"]
        )

    if (
        np.linalg.matrix_rank(gnsf["E"][np.ix_(ind_11, ind_11)]) != gnsf["nx1"]
        or z_check
    ):
        # print warning (always)
        print(f"the rank of E11 or E22 is not full after the reformulation")
        print("")
        print("the script will try to reformulate the model with an invertible matrix instead")
        print("NOTE: this feature is based on a heuristic, it should be used with care!!!")

        ## load models
        xdot = gnsf["xdot"]
        z = gnsf["z"]

        # # GNSF
        # get dimensions
        nx1 = gnsf["nx1"]
        x1dot = xdot[range(nx1)]

        k = ca.vertcat(x1dot, z)
        for i in [1, 2]:
            if i == 1:
                ind = range(gnsf["nx1"])
            else:
                ind = range(gnsf["nx1"], gnsf["nx1"] + gnsf["nz1"])
            mat = gnsf["E"][np.ix_(ind, ind)]
            breakpoint()

            while np.linalg.matrix_rank(mat) < len(ind):
                # import pdb; pdb.set_trace()
                if print_info:
                    print(" ")
                    print(f"the rank of E", str(i), str(i), " is not full")
                    print(
                        f"the algorithm will try to reformulate the model with an invertible matrix instead"
                    )
                    print(
                        f"NOTE: this feature is not super stable and might need more testing!!!!!!"
                    )
                for sub_max in ind:
                    sub_ind = range(min(ind), sub_max)
                    # regard the submatrix mat(sub_ind, sub_ind)
                    sub_mat = gnsf["E"][sub_ind, sub_ind]
                    if np.linalg.matrix_rank(sub_mat) < len(sub_ind):
                        # reformulate the model by adding a 1 to last diagonal
                        # element and changing rhs respectively.
                        gnsf["E"][sub_max, sub_max] = gnsf["E"][sub_max, sub_max] + 1
                        # this means adding the term 1 * k(sub_max) to the sub_max
                        # row of the l.h.s
                        if len(np.nonzero(gnsf["C"][sub_max, :])[0]) == 0:
                            # if isempty(find(gnsf['C'](sub_max,:), 1)):
                            # add new nonlinearity entry
                            gnsf["C"][sub_max, gnsf["n_out"] + 1] = 1
                            gnsf["phi_expr"] = ca.vertcat(gnsf["phi_expr"], k[sub_max])
                        else:
                            ind_f = np.nonzero(gnsf["C"][sub_max, :])[0]
                            if len(ind_f) != 1:
                                raise ValueError("C is assumed to be a selection matrix")
                            else:
                                ind_f = ind_f[0]
                            # add term to corresponding nonlinearity entry
                            # note: herbey we assume that C is a selection matrix,
                            # i.e. gnsf['phi_expr'](ind_f) is only entering one equation

                            gnsf["phi_expr"][ind_f] = (
                                gnsf["phi_expr"][ind_f]
                                + k[sub_max] / gnsf["C"][sub_max, ind_f]
                            )
                            gnsf = determine_input_nonlinearity_function(gnsf)
                            check_reformulation(model, gnsf, print_info)
        print("successfully reformulated the model with invertible matrices E11, E22")
    else:
        if print_info:
            print(" ")
            print(
                "the rank of both E11 and E22 is naturally full after the reformulation "
            )
            print("==>  model reformulation finished")
            print(" ")
    if (gnsf['nx2'] > 0 or gnsf['nz2'] > 0) and ca.det(gnsf["E_LO"]) == 0:
        print("_______________________________________")
        print(" ")
        print("TAKE CARE ")
        print("E_LO matrix is NOT regular after automatic transcription!")
        print("->> this means the model CANNOT be used with the gnsf integrator")
        print(
            "->> it probably means that one entry (of xdot or z) that was moved to the linear output type system"
        )
        print("    does not appear in the model at all (zero column in E_LO)")
        print(" OR: the columns of E_LO are linearly dependent ")
        print(" ")
        print(" SOLUTIONS: a) go through your model & check equations the method wanted to move to LOS")
        print("            b) deactivate the detect_LOS option")
        print("________________________________________")
    return gnsf



def reformulate_with_LOS(model: AcadosModel, gnsf, print_info):

    ## Description:
    # This function takes an intitial transcription of the implicit ODE model
    # "model" into "gnsf" and reformulates "gnsf" with a linear output system
    # (LOS), containing as many states of the model as possible.
    # Therefore it might be that the state vector and the implicit function
    # vector have to be reordered. This reordered model is part of the output,
    # namely reordered_model.

    # symbolics
    x = gnsf["x"]
    xdot = gnsf["xdot"]
    u = gnsf["u"]
    z = gnsf["z"]

    # dimensions
    nx = gnsf["nx"]
    nz = gnsf["nz"]

    # get model matrices
    A = gnsf["A"]
    B = gnsf["B"]
    C = gnsf["C"]
    E = gnsf["E"]
    c = gnsf["c"]

    A_LO = gnsf["A_LO"]

    y = gnsf["y"]

    phi_old = gnsf["phi_expr"]

    if print_info:
        print(" ")
        print("=================================================================")
        print("================    Detect Linear Output System   ===============")
        print("=================================================================")
        print(" ")
    ## build initial I_x1 and I_x2_candidates
    # I_xrange( all components of x for which either xii or xdot_ii enters y):
    # I_LOS_candidates: the remaining components

    I_nsf_components = set()
    I_LOS_candidates = set()

    if gnsf["ny"] > 0:
        for ii in range(nx):
            if ca.which_depends(y, x[ii])[0] or ca.which_depends(y, xdot[ii])[0]:
                # i.e. xii or xiidot are part of y, and enter phi_expr
                if print_info:
                    print(f"x_{ii} is part of x1")
                I_nsf_components = set.union(I_nsf_components, set([ii]))
            else:
                # i.e. neither xii nor xiidot are part of y, i.e. enter phi_expr
                I_LOS_candidates = set.union(I_LOS_candidates, set([ii]))
                if print_info:
                    print(" ")
        for ii in range(nz):
            if ca.which_depends(y, z[ii])[0]:
                # i.e. xii or xiidot are part of y, and enter phi_expr
                if print_info:
                    print(f"z_{ii} is part of x1")
                I_nsf_components = set.union(I_nsf_components, set([ii + nx]))
            else:
                # i.e. neither xii nor xiidot are part of y, i.e. enter phi_expr
                I_LOS_candidates = set.union(I_LOS_candidates, set([ii + nx]))
    else:
        I_LOS_candidates = set(range((nx + nz)))
    if print_info:
        print(" ")
        print(f"I_LOS_candidates {I_LOS_candidates}")
    new_nsf_components = I_nsf_components
    I_nsf_eq = set([])
    unsorted_dyn = set(range(nx + nz))
    xdot_z = ca.vertcat(xdot, z)

    ## determine components of Linear Output System
    # determine maximal index set I_x2
    # such that the components x(I_x2) can be written as a LOS
    Eq_map = []
    while True:
        ## find equations corresponding to new_nsf_components
        for ii in new_nsf_components:
            current_var = xdot_z[ii]
            var_name = current_var.name

            # print( unsorted_dyn)
            # print("np.nonzero(E[:,ii])[0]",np.nonzero(E[:,ii])[0])
            I_eq = set.intersection(set(np.nonzero(E[:, ii])[0]), unsorted_dyn)
            if len(I_eq) == 1:
                i_eq = I_eq.pop()
                if print_info:
                    print(f"component {i_eq} is associated with state {ii}")
            elif len(I_eq) > 1:  # x_ii_dot occurs in more than 1 eq linearly
                # find the equation with least linear dependencies on
                # I_LOS_cancidates
                number_of_eq = 0
                candidate_dependencies = np.zeros(len(I_eq), 1)
                I_x2_candidates = set.intersection(I_LOS_candidates, set(range(nx)))
                for eq in I_eq:
                    depending_candidates = set.union(
                        np.nonzero(E[eq, I_LOS_candidates])[0],
                        np.nonzero(A[eq, I_x2_candidates])[0],
                    )
                    candidate_dependencies[number_of_eq] = +len(depending_candidates)
                    number_of_eq += 1
                    number_of_eq = np.argmin(candidate_dependencies)
                i_eq = I_eq[number_of_eq]
            else:  ## x_ii_dot does not occur linearly in any of the unsorted dynamics
                for j in unsorted_dyn:
                    phi_eq_j = gnsf["phi_expr"][np.nonzero(C[j, :])[0]]
                    if ca.which_depends(phi_eq_j, xdot_z(ii))[0]:
                        I_eq = set.union(I_eq, j)
                if is_empty(I_eq):
                    I_eq = unsorted_dyn
                # find the equation with least linear dependencies on I_LOS_cancidates
                number_of_eq = 0
                candidate_dependencies = np.zeros(len(I_eq), 1)
                I_x2_candidates = set.intersection(I_LOS_candidates, set(range(nx)))
                for eq in I_eq:
                    depending_candidates = set.union(
                        np.nonzero(E[eq, I_LOS_candidates])[0],
                        np.nonzero(A[eq, I_x2_candidates])[0],
                    )
                    candidate_dependencies[number_of_eq] = +len(depending_candidates)
                    number_of_eq += 1
                    number_of_eq = np.argmin(candidate_dependencies)
                i_eq = I_eq[number_of_eq]
                ## add 1 * [xdot,z](ii) to both sides of i_eq
                if print_info:
                    print(f"adding 1 * {var_name} to both sides of equation {i_eq}.")
                gnsf["E"][i_eq, ii] = 1
                i_phi = np.nonzero(gnsf["C"][i_eq, :])
                if is_empty(i_phi):
                    i_phi = len(gnsf["phi_expr"]) + 1
                    gnsf["C"][i_eq, i_phi] = 1  # add column to C with 1 entry
                    gnsf["phi_expr"] = ca.vertcat(gnsf["phi_expr"], 0)
                    gnsf["phi_expr"][i_phi] = (
                        gnsf["phi_expr"](i_phi)
                        + gnsf["E"][i_eq, ii] / gnsf["C"][i_eq, i_phi] * xdot_z[ii]
                    )
                if print_info:
                    print(f"detected equation {i_eq} to correspond to variable {var_name}")
            I_nsf_eq = set.union(I_nsf_eq, {i_eq})
            # remove i_eq from unsorted_dyn
            unsorted_dyn.remove(i_eq)
            Eq_map.append([ii, i_eq])

        ## add components to I_x1
        for eq in I_nsf_eq:
            I_linear_dependence = set.union(
                set(np.nonzero(A[eq, :])[0]), set(np.nonzero(E[eq, :])[0])
            )
            I_nsf_components = set.union(I_linear_dependence, I_nsf_components)
            # I_nsf_components = I_nsf_components[:]

        new_nsf_components = set.intersection(I_LOS_candidates, I_nsf_components)
        if is_empty(new_nsf_components):
            if print_info:
                print("new_nsf_components is empty")
            break
        # remove new_nsf_components from candidates
        I_LOS_candidates = set.difference(I_LOS_candidates, new_nsf_components)
    if not is_empty(Eq_map):
        # [~, new_eq_order] = sort(Eq_map(1,:))
        # I_nsf_eq = Eq_map(2, new_eq_order )
        for count, m in enumerate(Eq_map):
            m.append(count)
        sorted(Eq_map, key=lambda x: x[1])
        new_eq_order = [m[2] for m in Eq_map]
        Eq_map = [Eq_map[i] for i in new_eq_order]
        I_nsf_eq = [m[1] for m in Eq_map]

    else:
        I_nsf_eq = []

    I_LOS_components = I_LOS_candidates
    I_LOS_eq = sorted(set.difference(set(range(nx + nz)), I_nsf_eq))
    I_nsf_eq = sorted(I_nsf_eq)

    I_x1 = set.intersection(I_nsf_components, set(range(nx)))
    I_z1 = set.intersection(I_nsf_components, set(range(nx, nx + nz)))
    I_z1 = set([i - nx for i in I_z1])

    I_x2 = set.intersection(I_LOS_components, set(range(nx)))
    I_z2 = set.intersection(I_LOS_components, set(range(nx, nx + nz)))
    I_z2 = set([i - nx for i in I_z2])

    if print_info:
        print(f"I_x1 {I_x1}, I_x2 {I_x2}")

    ## permute x, xdot
    if is_empty(I_x1):
        x1 = []
        x1dot = []
    else:
        x1 = x[list(I_x1)]
        x1dot = xdot[list(I_x1)]
    if is_empty(I_x2):
        x2 = []
        x2dot = []
    else:
        x2 = x[list(I_x2)]
        x2dot = xdot[list(I_x2)]
    if is_empty(I_z1):
        z1 = []
    else:
        z1 = z(I_z1)
    if is_empty(I_z2):
        z2 = []
    else:
        z2 = z[list(I_z2)]

    I_x1 = sorted(I_x1)
    I_x2 = sorted(I_x2)
    I_z1 = sorted(I_z1)
    I_z2 = sorted(I_z2)
    gnsf["xdot"] = ca.vertcat(x1dot, x2dot)
    gnsf["x"] = ca.vertcat(x1, x2)
    gnsf["z"] = ca.vertcat(z1, z2)
    gnsf["nx1"] = len(I_x1)
    gnsf["nx2"] = len(I_x2)
    gnsf["nz1"] = len(I_z1)
    gnsf["nz2"] = len(I_z2)

    # store permutations
    gnsf["idx_perm_x"] = I_x1 + I_x2
    gnsf["idx_perm_z"] = I_z1 + I_z2
    gnsf["idx_perm_f"] = I_nsf_eq + I_LOS_eq

    f_LO = ca.SX.sym("f_LO", 0, 0)

    ## rewrite I_LOS_eq as LOS
    if gnsf["n_out"] == 0:
        C_phi = np.zeros(gnsf["nx"] + gnsf["nz"], 1)
    else:
        C_phi = C @ phi_old
    if gnsf["nx1"] == 0:
        Ax1 = np.zeros(gnsf["nx"] + gnsf["nz"], 1)
    else:
        Ax1 = A[:, sorted(I_x1)] @ x1
    if gnsf["nx1"] + gnsf["nz1"] == 0:
        lhs_nsf = np.zeros(gnsf["nx"] + gnsf["nz"], 1)
    else:
        lhs_nsf = E[:, sorted(I_nsf_components)] @ ca.vertcat(x1, z1)
    n_LO = len(I_LOS_eq)
    B_LO = np.zeros((n_LO, gnsf["nu"]))
    A_LO = np.zeros((gnsf["nx2"] + gnsf["nz2"], gnsf["nx2"]))
    E_LO = np.zeros((n_LO, n_LO))
    c_LO = np.zeros((n_LO, 1))

    I_LOS_eq = list(I_LOS_eq)
    for eq in I_LOS_eq:
        i_LO = I_LOS_eq.index(eq)
        f_LO = ca.vertcat(f_LO, Ax1[eq] + C_phi[eq] - lhs_nsf[eq])
        print(f"eq {eq} I_LOS_components {I_LOS_components}, i_LO {i_LO}, f_LO {f_LO}")
        E_LO[i_LO, :] = E[eq, sorted(I_LOS_components)]
        A_LO[i_LO, :] = A[eq, I_x2]
        c_LO[i_LO, :] = c[eq]
        B_LO[i_LO, :] = B[eq, :]
    if casadi_length(f_LO) == 0:
        f_LO = ca.SX.zeros((gnsf["nx2"] + gnsf["nz2"], 1))
    f_LO = ca.simplify(f_LO)
    gnsf["A_LO"] = A_LO
    gnsf["E_LO"] = E_LO
    gnsf["B_LO"] = B_LO
    gnsf["c_LO"] = c_LO
    gnsf["f_lo_expr"] = f_LO

    ## remove I_LOS_eq from NSF type system
    gnsf["A"] = gnsf["A"][np.ix_(sorted(I_nsf_eq), sorted(I_x1))]
    gnsf["B"] = gnsf["B"][sorted(I_nsf_eq), :]
    gnsf["C"] = gnsf["C"][sorted(I_nsf_eq), :]
    gnsf["E"] = gnsf["E"][np.ix_(sorted(I_nsf_eq), sorted(I_nsf_components))]
    gnsf["c"] = gnsf["c"][sorted(I_nsf_eq), :]

    ## reduce phi, C
    I_nonzero = []
    for ii in range(gnsf["C"].shape[1]):  # n_colums of C:
        print(f"ii {ii}")
        if not all(gnsf["C"][:, ii] == 0):  # if column ~= 0
            I_nonzero.append(ii)
    gnsf["C"] = gnsf["C"][:, I_nonzero]
    gnsf["phi_expr"] = gnsf["phi_expr"][I_nonzero]

    gnsf = determine_input_nonlinearity_function(gnsf)

    check_reformulation(model, gnsf, print_info)

    gnsf["nontrivial_f_LO"] = 0
    if not is_empty(gnsf["f_lo_expr"]):
        for ii in range(casadi_length(gnsf["f_lo_expr"])):
            fii = gnsf["f_lo_expr"][ii]
            if not fii.is_zero():
                gnsf["nontrivial_f_LO"] = 1
            if not gnsf["nontrivial_f_LO"] and print_info:
                print("f_LO is fully trivial (== 0)")
    check_reformulation(model, gnsf, print_info)

    if print_info:
        print("")
        print("--------------------------------------------------------------")
        print("----- Success: Linear Output System (LOS) detected -----------")
        print("--------------------------------------------------------------")
        print("")
        print(f"==>>  moved  {gnsf['nx2']} differential states and {gnsf['nz2']} algebraic variables to the Linear Output System")
        print(f"==>>  reduced output dimension of phi from  {casadi_length(phi_old)} to {casadi_length(gnsf['phi_expr'])}")
        print(" ")
        print("Matrices defining the LOS read as")
        print(" ")
        print("E_LO =")
        print(gnsf["E_LO"])
        print("A_LO =")
        print(gnsf["A_LO"])
        print("B_LO =")
        print(gnsf["B_LO"])
        print("c_LO =")
        print(gnsf["c_LO"])

    return gnsf




def structure_detection_print_summary(gnsf, model):

    ## Description
    # this function prints the most important info after determining a GNSF
    # reformulation of the implicit model "initial_model" into "gnsf", which is
    # equivalent to the "reordered_model".
    # # GNSF
    # get dimensions
    nx = gnsf["nx"]
    nu = gnsf["nu"]
    nz = gnsf["nz"]

    nx1 = gnsf["nx1"]
    nx2 = gnsf["nx2"]

    nz1 = gnsf["nz1"]
    nz2 = gnsf["nz2"]

    # np = gnsf['np']
    n_out = gnsf["n_out"]
    ny = gnsf["ny"]
    nuhat = gnsf["nuhat"]

    n_nodes_initial = ca.n_nodes(model.f_impl_expr)
    x = gnsf["x"]
    z = gnsf["z"]

    phi_current = gnsf["phi_expr"]

    ## PRINT SUMMARY -- STRUCHTRE DETECTION
    print(" ")
    print("*************************************************************")
    print("****   SUCCESS: GNSF STRUCTURE DETECTION COMPLETE !!! *******")
    print("*************************************************************")
    print(" ")
    print(" ")
    print("=========== STRUCTURE DETECTION SUMMARY ======================"
    )
    print(" ")
    print("-------- Nonlinear Static Feedback type system --------")
    print(" ")
    print(f"reduced dimension of nonlinearity phi from        {nx + nz} to {gnsf['n_out']}")
    print(" ")
    print(f"reduced input dimension of nonlinearity phi from {2 * nx + nz + nu} to {gnsf['ny'] + gnsf['nuhat']}")
    print(" ")
    print(f"reduced number of nodes in CasADi expression of nonlinearity phi from  {n_nodes_initial}  to  {ca.n_nodes(phi_current)}\n")
    print("----------- Linear Output System (LOS) ---------------")
    if nx2 + nz2 > 0:
        print(" ")
        print(f"introduced Linear Output System of size           {nx2 + nz2}")
        print(" ")
        if nx2 > 0:
            print("consisting of the states:")
            print(" ")
            print(x[range(nx1, nx)])
            print(" ")
        if nz2 > 0:
            print("and algebraic variables:")
            print(" ")
            print(z[range(nz1, nz)])
            print(" ")
        if gnsf["purely_linear"] == 1:
            print(" ")
            print("Model is fully linear!")
            print(" ")
    if not all(gnsf["idx_perm_x"] == np.array(range(nx))):
        print(" ")
        print("----------------------------------------------------")
        print(
            "NOTE: permuted differential state vector x, such that x_gnsf = x(idx_perm_x) with idx_perm_x ="
        )
        print(" ")
        print(gnsf["idx_perm_x"])
    if nz != 0 and not all(gnsf["idx_perm_z"] == np.array(range(nz))):
        print(" ")
        print("----------------------------------------------------")
        print("NOTE: permuted algebraic state vector z, such that z_gnsf = z(idx_perm_z) with idx_perm_z =")
        print(" ")
        print(gnsf["idx_perm_z"])
    if not all(gnsf["idx_perm_f"] == np.array(range(nx + nz))):
        print(" ")
        print("----------------------------------------------------")
        print("NOTE: permuted rhs expression vector f, such that f_gnsf = f(idx_perm_f) with idx_perm_f =")
        print(" ")
        print(gnsf["idx_perm_f"])
    ## print GNSF dimensions
    print("----------------------------------------------------")
    print(" ")
    print("The dimensions of the GNSF reformulated model read as:")
    print(" ")
    print(f"nx    ", {nx})
    print(f"nu    ", {nu})
    print(f"nz    ", {nz})
    # print(f"np    ", {np})
    print(f"nx1   ", {nx1})
    print(f"nz1   ", {nz1})
    print(f"n_out ", {n_out})
    print(f"ny    ", {ny})
    print(f"nuhat ", {nuhat})


def idx_perm_to_ipiv(idx_perm):
    n = len(idx_perm)
    vec = list(range(n))
    ipiv = np.zeros(n)

    for ii in range(n):
        idx0 = idx_perm[ii]
        for jj in range(ii,n):
            if vec[jj]==idx0:
                idx1 = jj
                break
        tmp = vec[ii]
        vec[ii] = vec[idx1]
        vec[idx1] = tmp
        ipiv[ii] = idx1

    return ipiv

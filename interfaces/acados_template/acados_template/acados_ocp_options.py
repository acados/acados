# -*- coding: future_fstrings -*-
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

import numpy as np


INTEGRATOR_TYPES = ('ERK', 'IRK', 'GNSF', 'DISCRETE', 'LIFTED_IRK')
COLLOCATION_TYPES = ('GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE', 'EXPLICIT_RUNGE_KUTTA')
COST_DISCRETIZATION_TYPES = ('EULER', 'INTEGRATOR')


class AcadosOcpOptions:
    """
    class containing the description of the solver options
    """
    def __init__(self):
        self.__qp_solver = 'PARTIAL_CONDENSING_HPIPM'
        self.__hessian_approx = 'GAUSS_NEWTON'
        self.__integrator_type = 'ERK'
        self.__tf = None
        self.__nlp_solver_type = 'SQP_RTI'
        self.__globalization = 'FIXED_STEP'
        self.__nlp_solver_step_length = 1.0
        self.__levenberg_marquardt = 0.0
        self.__collocation_type = 'GAUSS_LEGENDRE'
        self.__sim_method_num_stages = 4
        self.__sim_method_num_steps = 1
        self.__sim_method_newton_iter = 3
        self.__sim_method_newton_tol = 0.0
        self.__sim_method_jac_reuse = 0
        self.__qp_solver_tol_stat = None
        self.__qp_solver_tol_eq = None
        self.__qp_solver_tol_ineq = None
        self.__qp_solver_tol_comp = None
        self.__qp_solver_iter_max = 50
        self.__qp_solver_cond_N = None
        self.__qp_solver_warm_start = 0
        self.__qp_solver_cond_ric_alg = 1
        self.__qp_solver_ric_alg = 1
        self.__nlp_solver_tol_stat = 1e-6
        self.__nlp_solver_tol_eq = 1e-6
        self.__nlp_solver_tol_ineq = 1e-6
        self.__nlp_solver_tol_comp = 1e-6
        self.__nlp_solver_max_iter = 100
        self.__nlp_solver_ext_qp_res = 0
        self.__Tsim = None
        self.__print_level = 0
        self.__initialize_t_slacks = 0
        self.__cost_discretization = 'EULER'
        self.__regularize_method = 'NO_REGULARIZE'
        self.__reg_epsilon = 1e-4
        self.__time_steps = None
        self.__shooting_nodes = None
        self.__exact_hess_cost = 1
        self.__exact_hess_dyn = 1
        self.__exact_hess_constr = 1
        self.__fixed_hess = 0
        self.__ext_cost_num_hess = 0
        self.__alpha_min = 0.05
        self.__alpha_reduction = 0.7
        self.__line_search_use_sufficient_descent = 0
        self.__globalization_use_SOC = 0
        self.__full_step_dual = 0
        self.__eps_sufficient_descent = 1e-4
        self.__hpipm_mode = 'BALANCE'
        # TODO: move those out? they are more about generation than about the acados OCP solver.
        self.__ext_fun_compile_flags = '-O2'
        self.__model_external_shared_lib_dir = None
        self.__model_external_shared_lib_name = None
        self.__custom_update_filename = ''
        self.__custom_update_header_filename = ''
        self.__custom_templates = []
        self.__custom_update_copy = True

    @property
    def qp_solver(self):
        """QP solver to be used in the NLP solver.
        String in ('PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP', 'FULL_CONDENSING_DAQP').
        Default: 'PARTIAL_CONDENSING_HPIPM'.
        """
        return self.__qp_solver

    @property
    def ext_fun_compile_flags(self):
        """
        String with compiler flags for external function compilation.
        Default: '-O2'.
        """
        return self.__ext_fun_compile_flags


    @property
    def custom_update_filename(self):
        """
        Filename of the custom C function to update solver data and parameters in between solver calls

        This file has to implement the functions
        int custom_update_init_function([model.name]_solver_capsule* capsule);
        int custom_update_function([model.name]_solver_capsule* capsule);
        int custom_update_terminate_function([model.name]_solver_capsule* capsule);


        Default: ''.
        """
        return self.__custom_update_filename


    @property
    def custom_templates(self):
        """
        List of tuples of the form:
        (input_filename, output_filename)

        Custom templates are render in OCP solver generation.

        Default: [].
        """
        return self.__custom_templates


    @property
    def custom_update_header_filename(self):
        """
        Header filename of the custom C function to update solver data and parameters in between solver calls

        This file has to declare the custom_update functions and look as follows:

        `// Called at the end of solver creation.`

        `// This is allowed to allocate memory and store the pointer to it into capsule->custom_update_memory.`

        `int custom_update_init_function([model.name]_solver_capsule* capsule);`

        `// Custom update function that can be called between solver calls`

        `int custom_update_function([model.name]_solver_capsule* capsule, double* data, int data_len);`

        `// Called just before destroying the solver.`

        `// Responsible to free allocated memory, stored at capsule->custom_update_memory.`

        `int custom_update_terminate_function([model.name]_solver_capsule* capsule);`

        Default: ''.
        """
        return self.__custom_update_header_filename

    @property
    def custom_update_copy(self):
        """
        Boolean;
        If True, the custom update function files are copied into the `code_export_directory`.
        """
        return self.__custom_update_copy


    @property
    def hpipm_mode(self):
        """
        Mode of HPIPM to be used,

        String in ('BALANCE', 'SPEED_ABS', 'SPEED', 'ROBUST').

        Default: 'BALANCE'.

        see https://cdn.syscop.de/publications/Frison2020a.pdf
        and the HPIPM code:
        https://github.com/giaf/hpipm/blob/master/ocp_qp/x_ocp_qp_ipm.c#L69
        """
        return self.__hpipm_mode

    @property
    def hessian_approx(self):
        """Hessian approximation.
        String in ('GAUSS_NEWTON', 'EXACT').
        Default: 'GAUSS_NEWTON'.
        """
        return self.__hessian_approx

    @property
    def integrator_type(self):
        """
        Integrator type.
        String in ('ERK', 'IRK', 'GNSF', 'DISCRETE', 'LIFTED_IRK').
        Default: 'ERK'.
        """
        return self.__integrator_type

    @property
    def nlp_solver_type(self):
        """NLP solver.
        String in ('SQP', 'SQP_RTI').
        Default: 'SQP_RTI'.
        """
        return self.__nlp_solver_type

    @property
    def globalization(self):
        """Globalization type.
        String in ('FIXED_STEP', 'MERIT_BACKTRACKING').
        Default: 'FIXED_STEP'.

        .. note:: preliminary implementation.
        """
        return self.__globalization

    @property
    def collocation_type(self):
        """Collocation type: only relevant for implicit integrators
        -- string in {'GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE', 'EXPLICIT_RUNGE_KUTTA'}.

        Default: GAUSS_LEGENDRE
        """
        return self.__collocation_type

    @property
    def regularize_method(self):
        """Regularization method for the Hessian.
        String in ('NO_REGULARIZE', 'MIRROR', 'PROJECT', 'PROJECT_REDUC_HESS', 'CONVEXIFY') or :code:`None`.

        - MIRROR: performs eigenvalue decomposition H = V^T D V and sets D_ii = max(eps, abs(D_ii))
        - PROJECT: performs eigenvalue decomposition H = V^T D V and sets D_ii = max(eps, D_ii)
        - CONVEXIFY: Algorithm 6 from Verschueren2017, https://cdn.syscop.de/publications/Verschueren2017.pdf
        - PROJECT_REDUC_HESS: experimental

        Note: default eps = 1e-4

        Default: 'NO_REGULARIZE'.
        """
        return self.__regularize_method

    @property
    def nlp_solver_step_length(self):
        """
        Fixed Newton step length.
        Type: float >= 0.
        Default: 1.0.
        """
        return self.__nlp_solver_step_length

    @property
    def levenberg_marquardt(self):
        """
        Factor for LM regularization.
        Type: float >= 0
        Default: 0.0.
        """
        return self.__levenberg_marquardt

    @property
    def sim_method_num_stages(self):
        """
        Number of stages in the integrator.
        Type: int > 0 or ndarray of ints > 0 of shape (N,).
        Default: 4
        """
        return self.__sim_method_num_stages

    @property
    def sim_method_num_steps(self):
        """
        Number of steps in the integrator.
        Type: int > 0 or ndarray of ints > 0 of shape (N,).
        Default: 1
        """
        return self.__sim_method_num_steps

    @property
    def sim_method_newton_iter(self):
        """
        Number of Newton iterations in simulation method.
        Type: int > 0
        Default: 3
        """
        return self.__sim_method_newton_iter

    @property
    def sim_method_newton_tol(self):
        """
        Tolerance of Newton system in simulation method.
        Type: float: 0.0 means not used
        Default: 0.0
        """
        return self.__sim_method_newton_tol

    @property
    def sim_method_jac_reuse(self):
        """
        Integer determining if jacobians are reused within integrator or ndarray of ints > 0 of shape (N,).
        0: False (no reuse); 1: True (reuse)
        Default: 0
        """
        return self.__sim_method_jac_reuse

    @property
    def qp_solver_tol_stat(self):
        """
        QP solver stationarity tolerance.
        Default: :code:`None`
        """
        return self.__qp_solver_tol_stat

    @property
    def qp_solver_tol_eq(self):
        """
        QP solver equality tolerance.
        Default: :code:`None`
        """
        return self.__qp_solver_tol_eq

    @property
    def qp_solver_tol_ineq(self):
        """
        QP solver inequality.
        Default: :code:`None`
        """
        return self.__qp_solver_tol_ineq

    @property
    def qp_solver_tol_comp(self):
        """
        QP solver complementarity.
        Default: :code:`None`
        """
        return self.__qp_solver_tol_comp

    @property
    def qp_solver_cond_N(self):
        """QP solver: New horizon after partial condensing.
        Set to N by default -> no condensing."""
        return self.__qp_solver_cond_N

    @property
    def qp_solver_warm_start(self):
        """
        QP solver: Warm starting.
        0: no warm start; 1: warm start; 2: hot start.
        Default: 0
        """
        return self.__qp_solver_warm_start

    @property
    def qp_solver_cond_ric_alg(self):
        """
        QP solver: Determines which algorithm is used in HPIPM condensing.
        0: dont factorize hessian in the condensing; 1: factorize.
        Default: 1
        """
        return self.__qp_solver_cond_ric_alg

    @property
    def qp_solver_ric_alg(self):
        """
        QP solver: Determines which algorithm is used in HPIPM OCP QP solver.
        0 classical Riccati, 1 square-root Riccati.

        Note: taken from [HPIPM paper]:

        (a) the classical implementation requires the reduced  (projected) Hessian with respect to the dynamics
            equality constraints to be positive definite, but allows the full-space Hessian to be indefinite)
        (b) the square-root implementation, which in order to reduce the flop count employs the Cholesky
            factorization of the Riccati recursion matrix (P), and therefore requires the full-space Hessian to be positive definite

        [HPIPM paper]: HPIPM: a high-performance quadratic programming framework for model predictive control, Frison and Diehl, 2020
        https://cdn.syscop.de/publications/Frison2020a.pdf

        Default: 1
        """
        return self.__qp_solver_ric_alg

    @property
    def qp_solver_iter_max(self):
        """
        QP solver: maximum number of iterations.
        Type: int > 0
        Default: 50
        """
        return self.__qp_solver_iter_max

    @property
    def tol(self):
        """
        NLP solver tolerance. Sets or gets the max of :py:attr:`nlp_solver_tol_eq`,
        :py:attr:`nlp_solver_tol_ineq`, :py:attr:`nlp_solver_tol_comp`
        and :py:attr:`nlp_solver_tol_stat`.
        """
        return max([self.__nlp_solver_tol_eq, self.__nlp_solver_tol_ineq,\
                    self.__nlp_solver_tol_comp, self.__nlp_solver_tol_stat])

    @property
    def qp_tol(self):
        """
        QP solver tolerance.
        Sets all of the following at once or gets the max of
        :py:attr:`qp_solver_tol_eq`, :py:attr:`qp_solver_tol_ineq`,
        :py:attr:`qp_solver_tol_comp` and
        :py:attr:`qp_solver_tol_stat`.
        """
        return max([self.__qp_solver_tol_eq, self.__qp_solver_tol_ineq,\
                    self.__qp_solver_tol_comp, self.__qp_solver_tol_stat])

    @property
    def nlp_solver_tol_stat(self):
        """
        NLP solver stationarity tolerance.
        Type: float > 0
        Default: 1e-6
        """
        return self.__nlp_solver_tol_stat

    @property
    def nlp_solver_tol_eq(self):
        """NLP solver equality tolerance"""
        return self.__nlp_solver_tol_eq

    @property
    def alpha_min(self):
        """Minimal step size for globalization MERIT_BACKTRACKING, default: 0.05."""
        return self.__alpha_min

    @property
    def reg_epsilon(self):
        """Epsilon for regularization, used if regularize_method in ['PROJECT', 'MIRROR', 'CONVEXIFY']"""
        return self.__reg_epsilon

    @property
    def alpha_reduction(self):
        """Step size reduction factor for globalization MERIT_BACKTRACKING, default: 0.7."""
        return self.__alpha_reduction

    @property
    def line_search_use_sufficient_descent(self):
        """
        Determines if sufficient descent (Armijo) condition is used in line search.
        Type: int; 0 or 1;
        default: 0.
        """
        return self.__line_search_use_sufficient_descent

    @property
    def eps_sufficient_descent(self):
        """
        Factor for sufficient descent (Armijo) conditon, see line_search_use_sufficient_descent.
        Type: float,
        default: 1e-4.
        """
        return self.__eps_sufficient_descent

    @property
    def globalization_use_SOC(self):
        """
        Determines if second order correction (SOC) is done when using MERIT_BACKTRACKING.
        SOC is done if preliminary line search does not return full step.
        Type: int; 0 or 1;
        default: 0.
        """
        return self.__globalization_use_SOC

    @property
    def full_step_dual(self):
        """
        Determines if dual variables are updated with full steps (alpha=1.0) when primal variables are updated with smaller step.
        Type: int; 0 or 1;
        default: 0.
        """
        return self.__full_step_dual

    @property
    def nlp_solver_tol_ineq(self):
        """NLP solver inequality tolerance"""
        return self.__nlp_solver_tol_ineq

    @property
    def nlp_solver_ext_qp_res(self):
        """Determines if residuals of QP are computed externally within NLP solver (for debugging)

        Type: int; 0 or 1;
        Default: 0.
        """
        return self.__nlp_solver_ext_qp_res

    @property
    def nlp_solver_tol_comp(self):
        """NLP solver complementarity tolerance"""
        return self.__nlp_solver_tol_comp

    @property
    def nlp_solver_max_iter(self):
        """
        NLP solver maximum number of iterations.
        Type: int > 0
        Default: 100
        """
        return self.__nlp_solver_max_iter

    @property
    def time_steps(self):
        """
        Vector with time steps between the shooting nodes. Set automatically to uniform discretization if :py:attr:`N` and :py:attr:`tf` are provided.
        Default: :code:`None`
        """
        return self.__time_steps

    @property
    def shooting_nodes(self):
        """
        Vector with the shooting nodes, time_steps will be computed from it automatically.
        Default: :code:`None`
        """
        return self.__shooting_nodes

    @property
    def tf(self):
        """
        Prediction horizon
        Type: float > 0
        Default: :code:`None`
        """
        return self.__tf

    @property
    def Tsim(self):
        """
        Time horizon for one integrator step. Automatically calculated as :py:attr:`tf`/:py:attr:`N`.
        Default: :code:`None`
        """
        return self.__Tsim

    @property
    def print_level(self):
        """
        Verbosity of printing.
        Type: int >= 0
        Default: 0
        """
        return self.__print_level

    @property
    def model_external_shared_lib_dir(self):
        """Path to the .so lib"""
        return self.__model_external_shared_lib_dir

    @property
    def model_external_shared_lib_name(self):
        """Name of the .so lib"""
        return self.__model_external_shared_lib_name

    @property
    def exact_hess_constr(self):
        """
        Used in case of hessian_approx == 'EXACT'.\n
        Can be used to turn off exact hessian contributions from the constraints module.
        """
        return self.__exact_hess_constr

    @property
    def exact_hess_cost(self):
        """
        Used in case of hessian_approx == 'EXACT'.\n
        Can be used to turn off exact hessian contributions from the cost module.
        """
        return self.__exact_hess_cost

    @property
    def exact_hess_dyn(self):
        """
        Used in case of hessian_approx == 'EXACT'.\n
        Can be used to turn off exact hessian contributions from the dynamics module.
        """
        return self.__exact_hess_dyn

    @property
    def fixed_hess(self):
        """
        Indicates if the hessian is fixed (1) or not (0).\n
        If fixed, the hessian is computed only once and not updated in the SQP loop.
        This can safely be set to 1 if there are no slacked constraints, cost module is 'LINEAR_LS' with Gauss-Newton Hessian,
        and the weighting matrix 'W' is not updated after solver creation.
        Default: 0
        """
        return self.__fixed_hess

    @property
    def ext_cost_num_hess(self):
        """
        Determines if custom hessian approximation for cost contribution is used (> 0).\n
        Or if hessian contribution is evaluated exactly using CasADi external function (=0 - default).
        """
        return self.__ext_cost_num_hess

    @property
    def cost_discretization(self):
        """
        Cost discretization: string in {'EULER', 'INTEGRATOR'}.
        Default: 'EULER'
        'EULER': cost is evaluated at shooting nodes
        'INTEGRATOR': cost is integrated over the shooting intervals - only supported for IRK integrator
        """
        return self.__cost_discretization

    @qp_solver.setter
    def qp_solver(self, qp_solver):
        qp_solvers = ('PARTIAL_CONDENSING_HPIPM', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM', \
                'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP', \
                'FULL_CONDENSING_DAQP')
        if qp_solver in qp_solvers:
            self.__qp_solver = qp_solver
        else:
            raise Exception('Invalid qp_solver value. Possible values are:\n\n' \
                    + ',\n'.join(qp_solvers) + '.\n\nYou have: ' + qp_solver + '.\n\n')

    @regularize_method.setter
    def regularize_method(self, regularize_method):
        regularize_methods = ('NO_REGULARIZE', 'MIRROR', 'PROJECT', \
                                'PROJECT_REDUC_HESS', 'CONVEXIFY')
        if regularize_method in regularize_methods:
            self.__regularize_method = regularize_method
        else:
            raise Exception('Invalid regularize_method value. Possible values are:\n\n' \
                    + ',\n'.join(regularize_methods) + '.\n\nYou have: ' + regularize_method + '.\n\n')

    @collocation_type.setter
    def collocation_type(self, collocation_type):
        if collocation_type in COLLOCATION_TYPES:
            self.__collocation_type = collocation_type
        else:
            raise Exception('Invalid collocation_type value. Possible values are:\n\n' \
                    + ',\n'.join(COLLOCATION_TYPES) + '.\n\nYou have: ' + collocation_type + '.\n\n')

    @hpipm_mode.setter
    def hpipm_mode(self, hpipm_mode):
        hpipm_modes = ('BALANCE', 'SPEED_ABS', 'SPEED', 'ROBUST')
        if hpipm_mode in hpipm_modes:
            self.__hpipm_mode = hpipm_mode
        else:
            raise Exception('Invalid hpipm_mode value. Possible values are:\n\n' \
                    + ',\n'.join(hpipm_modes) + '.\n\nYou have: ' + hpipm_mode + '.\n\n')

    @ext_fun_compile_flags.setter
    def ext_fun_compile_flags(self, ext_fun_compile_flags):
        if isinstance(ext_fun_compile_flags, str):
            self.__ext_fun_compile_flags = ext_fun_compile_flags
        else:
            raise Exception('Invalid ext_fun_compile_flags, expected a string.\n')


    @custom_update_filename.setter
    def custom_update_filename(self, custom_update_filename):
        if isinstance(custom_update_filename, str):
            self.__custom_update_filename = custom_update_filename
        else:
            raise Exception('Invalid custom_update_filename, expected a string.\n')

    @custom_templates.setter
    def custom_templates(self, custom_templates):
        if not isinstance(custom_templates, list):
            raise Exception('Invalid custom_templates, expected a list.\n')
        for tup in custom_templates:
            if not isinstance(tup, tuple):
                raise Exception('Invalid custom_templates, shoubld be list of tuples.\n')
            for s in tup:
                if not isinstance(s, str):
                    raise Exception('Invalid custom_templates, shoubld be list of tuples of strings.\n')
        self.__custom_templates = custom_templates

    @custom_update_header_filename.setter
    def custom_update_header_filename(self, custom_update_header_filename):
        if isinstance(custom_update_header_filename, str):
            self.__custom_update_header_filename = custom_update_header_filename
        else:
            raise Exception('Invalid custom_update_header_filename, expected a string.\n')

    @custom_update_copy.setter
    def custom_update_copy(self, custom_update_copy):
        if isinstance(custom_update_copy, bool):
            self.__custom_update_copy = custom_update_copy
        else:
            raise Exception('Invalid custom_update_copy, expected a bool.\n')

    @hessian_approx.setter
    def hessian_approx(self, hessian_approx):
        hessian_approxs = ('GAUSS_NEWTON', 'EXACT')
        if hessian_approx in hessian_approxs:
            self.__hessian_approx = hessian_approx
        else:
            raise Exception('Invalid hessian_approx value. Possible values are:\n\n' \
                    + ',\n'.join(hessian_approxs) + '.\n\nYou have: ' + hessian_approx + '.\n\n')

    @integrator_type.setter
    def integrator_type(self, integrator_type):
        if integrator_type in INTEGRATOR_TYPES:
            self.__integrator_type = integrator_type
        else:
            raise Exception('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(INTEGRATOR_TYPES) + '.\n\nYou have: ' + integrator_type + '.\n\n')

    @tf.setter
    def tf(self, tf):
        self.__tf = tf

    @time_steps.setter
    def time_steps(self, time_steps):
        if isinstance(time_steps, np.ndarray):
            if len(time_steps.shape) == 1:
                    self.__time_steps = time_steps
            else:
                raise Exception('Invalid time_steps, expected np.ndarray of shape (N,).')
        else:
            raise Exception('Invalid time_steps, expected np.ndarray.')

    @shooting_nodes.setter
    def shooting_nodes(self, shooting_nodes):
        if isinstance(shooting_nodes, np.ndarray):
            if len(shooting_nodes.shape) == 1:
                self.__shooting_nodes = shooting_nodes
            else:
                raise Exception('Invalid shooting_nodes, expected np.ndarray of shape (N+1,).')
        else:
            raise Exception('Invalid shooting_nodes, expected np.ndarray.')

    @Tsim.setter
    def Tsim(self, Tsim):
        self.__Tsim = Tsim

    @globalization.setter
    def globalization(self, globalization):
        globalization_types = ('MERIT_BACKTRACKING', 'FIXED_STEP')
        if globalization in globalization_types:
            self.__globalization = globalization
        else:
            raise Exception('Invalid globalization value. Possible values are:\n\n' \
                    + ',\n'.join(globalization_types) + '.\n\nYou have: ' + globalization + '.\n\n')

    @reg_epsilon.setter
    def reg_epsilon(self, reg_epsilon):
        self.__reg_epsilon = reg_epsilon

    @alpha_min.setter
    def alpha_min(self, alpha_min):
        self.__alpha_min = alpha_min

    @alpha_reduction.setter
    def alpha_reduction(self, alpha_reduction):
        self.__alpha_reduction = alpha_reduction

    @line_search_use_sufficient_descent.setter
    def line_search_use_sufficient_descent(self, line_search_use_sufficient_descent):
        if line_search_use_sufficient_descent in [0, 1]:
            self.__line_search_use_sufficient_descent = line_search_use_sufficient_descent
        else:
            raise Exception(f'Invalid value for line_search_use_sufficient_descent. Possible values are 0, 1, got {line_search_use_sufficient_descent}')

    @globalization_use_SOC.setter
    def globalization_use_SOC(self, globalization_use_SOC):
        if globalization_use_SOC in [0, 1]:
            self.__globalization_use_SOC = globalization_use_SOC
        else:
            raise Exception(f'Invalid value for globalization_use_SOC. Possible values are 0, 1, got {globalization_use_SOC}')

    @full_step_dual.setter
    def full_step_dual(self, full_step_dual):
        if full_step_dual in [0, 1]:
            self.__full_step_dual = full_step_dual
        else:
            raise Exception(f'Invalid value for full_step_dual. Possible values are 0, 1, got {full_step_dual}')

    @eps_sufficient_descent.setter
    def eps_sufficient_descent(self, eps_sufficient_descent):
        if isinstance(eps_sufficient_descent, float) and eps_sufficient_descent > 0:
            self.__eps_sufficient_descent = eps_sufficient_descent
        else:
            raise Exception('Invalid eps_sufficient_descent value. eps_sufficient_descent must be a positive float.')

    @sim_method_num_stages.setter
    def sim_method_num_stages(self, sim_method_num_stages):

        # if isinstance(sim_method_num_stages, int):
        #     self.__sim_method_num_stages = sim_method_num_stages
        # else:
        #     raise Exception('Invalid sim_method_num_stages value. sim_method_num_stages must be an integer.')

        self.__sim_method_num_stages = sim_method_num_stages

    @sim_method_num_steps.setter
    def sim_method_num_steps(self, sim_method_num_steps):

        # if isinstance(sim_method_num_steps, int):
        #     self.__sim_method_num_steps = sim_method_num_steps
        # else:
        #     raise Exception('Invalid sim_method_num_steps value. sim_method_num_steps must be an integer.')
        self.__sim_method_num_steps = sim_method_num_steps


    @sim_method_newton_iter.setter
    def sim_method_newton_iter(self, sim_method_newton_iter):

        if isinstance(sim_method_newton_iter, int):
            self.__sim_method_newton_iter = sim_method_newton_iter
        else:
            raise Exception('Invalid sim_method_newton_iter value. sim_method_newton_iter must be an integer.')

    @sim_method_newton_tol.setter
    def sim_method_newton_tol(self, sim_method_newton_tol):
        if isinstance(sim_method_newton_tol, float) and sim_method_newton_tol > 0:
            self.__sim_method_newton_tol = sim_method_newton_tol
        else:
            raise Exception('Invalid sim_method_newton_tol value. sim_method_newton_tol must be a positive float.')

    @sim_method_jac_reuse.setter
    def sim_method_jac_reuse(self, sim_method_jac_reuse):
        self.__sim_method_jac_reuse = sim_method_jac_reuse

    @nlp_solver_type.setter
    def nlp_solver_type(self, nlp_solver_type):
        nlp_solver_types = ('SQP', 'SQP_RTI')
        if nlp_solver_type in nlp_solver_types:
            self.__nlp_solver_type = nlp_solver_type
        else:
            raise Exception('Invalid nlp_solver_type value. Possible values are:\n\n' \
                    + ',\n'.join(nlp_solver_types) + '.\n\nYou have: ' + nlp_solver_type + '.\n\n')

    @cost_discretization.setter
    def cost_discretization(self, cost_discretization):
        if cost_discretization in COST_DISCRETIZATION_TYPES:
            self.__cost_discretization = cost_discretization
        else:
            raise Exception('Invalid cost_discretization value. Possible values are:\n\n' \
                    + ',\n'.join(COST_DISCRETIZATION_TYPES) + '.\n\nYou have: ' + cost_discretization + '.')

    @nlp_solver_step_length.setter
    def nlp_solver_step_length(self, nlp_solver_step_length):
        if isinstance(nlp_solver_step_length, float) and nlp_solver_step_length >= 0:
            self.__nlp_solver_step_length = nlp_solver_step_length
        else:
            raise Exception('Invalid nlp_solver_step_length value. nlp_solver_step_length must be a positive float.')

    @levenberg_marquardt.setter
    def levenberg_marquardt(self, levenberg_marquardt):
        if isinstance(levenberg_marquardt, float) and levenberg_marquardt >= 0:
            self.__levenberg_marquardt = levenberg_marquardt
        else:
            raise Exception('Invalid levenberg_marquardt value. levenberg_marquardt must be a positive float.')

    @qp_solver_iter_max.setter
    def qp_solver_iter_max(self, qp_solver_iter_max):
        if isinstance(qp_solver_iter_max, int) and qp_solver_iter_max > 0:
            self.__qp_solver_iter_max = qp_solver_iter_max
        else:
            raise Exception('Invalid qp_solver_iter_max value. qp_solver_iter_max must be a positive int.')

    @qp_solver_ric_alg.setter
    def qp_solver_ric_alg(self, qp_solver_ric_alg):
        if qp_solver_ric_alg in [0, 1]:
            self.__qp_solver_ric_alg = qp_solver_ric_alg
        else:
            raise Exception(f'Invalid qp_solver_ric_alg value. qp_solver_ric_alg must be in [0, 1], got {qp_solver_ric_alg}.')

    @qp_solver_cond_ric_alg.setter
    def qp_solver_cond_ric_alg(self, qp_solver_cond_ric_alg):
        if qp_solver_cond_ric_alg in [0, 1]:
            self.__qp_solver_cond_ric_alg = qp_solver_cond_ric_alg
        else:
            raise Exception(f'Invalid qp_solver_cond_ric_alg value. qp_solver_cond_ric_alg must be in [0, 1], got {qp_solver_cond_ric_alg}.')


    @qp_solver_cond_N.setter
    def qp_solver_cond_N(self, qp_solver_cond_N):
        if isinstance(qp_solver_cond_N, int) and qp_solver_cond_N >= 0:
            self.__qp_solver_cond_N = qp_solver_cond_N
        else:
            raise Exception('Invalid qp_solver_cond_N value. qp_solver_cond_N must be a positive int.')

    @qp_solver_warm_start.setter
    def qp_solver_warm_start(self, qp_solver_warm_start):
        if qp_solver_warm_start in [0, 1, 2]:
            self.__qp_solver_warm_start = qp_solver_warm_start
        else:
            raise Exception('Invalid qp_solver_warm_start value. qp_solver_warm_start must be 0 or 1 or 2.')

    @qp_tol.setter
    def qp_tol(self, qp_tol):
        if isinstance(qp_tol, float) and qp_tol > 0:
            self.__qp_solver_tol_eq = qp_tol
            self.__qp_solver_tol_ineq = qp_tol
            self.__qp_solver_tol_stat = qp_tol
            self.__qp_solver_tol_comp = qp_tol
        else:
            raise Exception('Invalid qp_tol value. qp_tol must be a positive float.')

    @qp_solver_tol_stat.setter
    def qp_solver_tol_stat(self, qp_solver_tol_stat):
        if isinstance(qp_solver_tol_stat, float) and qp_solver_tol_stat > 0:
            self.__qp_solver_tol_stat = qp_solver_tol_stat
        else:
            raise Exception('Invalid qp_solver_tol_stat value. qp_solver_tol_stat must be a positive float.')

    @qp_solver_tol_eq.setter
    def qp_solver_tol_eq(self, qp_solver_tol_eq):
        if isinstance(qp_solver_tol_eq, float) and qp_solver_tol_eq > 0:
            self.__qp_solver_tol_eq = qp_solver_tol_eq
        else:
            raise Exception('Invalid qp_solver_tol_eq value. qp_solver_tol_eq must be a positive float.')

    @qp_solver_tol_ineq.setter
    def qp_solver_tol_ineq(self, qp_solver_tol_ineq):
        if isinstance(qp_solver_tol_ineq, float) and qp_solver_tol_ineq > 0:
            self.__qp_solver_tol_ineq = qp_solver_tol_ineq
        else:
            raise Exception('Invalid qp_solver_tol_ineq value. qp_solver_tol_ineq must be a positive float.')

    @qp_solver_tol_comp.setter
    def qp_solver_tol_comp(self, qp_solver_tol_comp):
        if isinstance(qp_solver_tol_comp, float) and qp_solver_tol_comp > 0:
            self.__qp_solver_tol_comp = qp_solver_tol_comp
        else:
            raise Exception('Invalid qp_solver_tol_comp value. qp_solver_tol_comp must be a positive float.')

    @tol.setter
    def tol(self, tol):
        if isinstance(tol, float) and tol > 0:
            self.__nlp_solver_tol_eq = tol
            self.__nlp_solver_tol_ineq = tol
            self.__nlp_solver_tol_stat = tol
            self.__nlp_solver_tol_comp = tol
        else:
            raise Exception('Invalid tol value. tol must be a positive float.')

    @nlp_solver_tol_stat.setter
    def nlp_solver_tol_stat(self, nlp_solver_tol_stat):
        if isinstance(nlp_solver_tol_stat, float) and nlp_solver_tol_stat > 0:
            self.__nlp_solver_tol_stat = nlp_solver_tol_stat
        else:
            raise Exception('Invalid nlp_solver_tol_stat value. nlp_solver_tol_stat must be a positive float.')

    @nlp_solver_tol_eq.setter
    def nlp_solver_tol_eq(self, nlp_solver_tol_eq):
        if isinstance(nlp_solver_tol_eq, float) and nlp_solver_tol_eq > 0:
            self.__nlp_solver_tol_eq = nlp_solver_tol_eq
        else:
            raise Exception('Invalid nlp_solver_tol_eq value. nlp_solver_tol_eq must be a positive float.')

    @nlp_solver_tol_ineq.setter
    def nlp_solver_tol_ineq(self, nlp_solver_tol_ineq):
        if isinstance(nlp_solver_tol_ineq, float) and nlp_solver_tol_ineq > 0:
            self.__nlp_solver_tol_ineq = nlp_solver_tol_ineq
        else:
            raise Exception('Invalid nlp_solver_tol_ineq value. nlp_solver_tol_ineq must be a positive float.')

    @nlp_solver_ext_qp_res.setter
    def nlp_solver_ext_qp_res(self, nlp_solver_ext_qp_res):
        if nlp_solver_ext_qp_res in [0, 1]:
            self.__nlp_solver_ext_qp_res = nlp_solver_ext_qp_res
        else:
            raise Exception('Invalid nlp_solver_ext_qp_res value. nlp_solver_ext_qp_res must be in [0, 1].')

    @nlp_solver_tol_comp.setter
    def nlp_solver_tol_comp(self, nlp_solver_tol_comp):
        if isinstance(nlp_solver_tol_comp, float) and nlp_solver_tol_comp > 0:
            self.__nlp_solver_tol_comp = nlp_solver_tol_comp
        else:
            raise Exception('Invalid nlp_solver_tol_comp value. nlp_solver_tol_comp must be a positive float.')

    @nlp_solver_max_iter.setter
    def nlp_solver_max_iter(self, nlp_solver_max_iter):

        if isinstance(nlp_solver_max_iter, int) and nlp_solver_max_iter > 0:
            self.__nlp_solver_max_iter = nlp_solver_max_iter
        else:
            raise Exception('Invalid nlp_solver_max_iter value. nlp_solver_max_iter must be a positive int.')

    @print_level.setter
    def print_level(self, print_level):
        if isinstance(print_level, int) and print_level >= 0:
            self.__print_level = print_level
        else:
            raise Exception('Invalid print_level value. print_level takes one of the values >=0.')

    @model_external_shared_lib_dir.setter
    def model_external_shared_lib_dir(self, model_external_shared_lib_dir):
        if isinstance(model_external_shared_lib_dir, str) :
            self.__model_external_shared_lib_dir = model_external_shared_lib_dir
        else:
            raise Exception('Invalid model_external_shared_lib_dir value. Str expected.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_dir) + '.\n\n')

    @model_external_shared_lib_name.setter
    def model_external_shared_lib_name(self, model_external_shared_lib_name):
        if isinstance(model_external_shared_lib_name, str) :
            if model_external_shared_lib_name[-3:] == '.so' :
                raise Exception('Invalid model_external_shared_lib_name value. Remove the .so extension.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\n')
            else :
                self.__model_external_shared_lib_name = model_external_shared_lib_name
        else:
            raise Exception('Invalid model_external_shared_lib_name value. Str expected.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\n')

    @exact_hess_constr.setter
    def exact_hess_constr(self, exact_hess_constr):
        if exact_hess_constr in [0, 1]:
            self.__exact_hess_constr = exact_hess_constr
        else:
            raise Exception('Invalid exact_hess_constr value. exact_hess_constr takes one of the values 0, 1.')

    @exact_hess_cost.setter
    def exact_hess_cost(self, exact_hess_cost):
        if exact_hess_cost in [0, 1]:
            self.__exact_hess_cost = exact_hess_cost
        else:
            raise Exception('Invalid exact_hess_cost value. exact_hess_cost takes one of the values 0, 1.')

    @exact_hess_dyn.setter
    def exact_hess_dyn(self, exact_hess_dyn):
        if exact_hess_dyn in [0, 1]:
            self.__exact_hess_dyn = exact_hess_dyn
        else:
            raise Exception('Invalid exact_hess_dyn value. exact_hess_dyn takes one of the values 0, 1.')

    @fixed_hess.setter
    def fixed_hess(self, fixed_hess):
        if fixed_hess in [0, 1]:
            self.__fixed_hess = fixed_hess
        else:
            raise Exception('Invalid fixed_hess value. fixed_hess takes one of the values 0, 1.')

    @ext_cost_num_hess.setter
    def ext_cost_num_hess(self, ext_cost_num_hess):
        if ext_cost_num_hess in [0, 1]:
            self.__ext_cost_num_hess = ext_cost_num_hess
        else:
            raise Exception('Invalid ext_cost_num_hess value. ext_cost_num_hess takes one of the values 0, 1.')

    def set(self, attr, value):
        setattr(self, attr, value)

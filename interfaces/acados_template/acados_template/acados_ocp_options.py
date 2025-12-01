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

import os
import warnings

from deprecated.sphinx import deprecated
from .utils import check_if_nparray_and_flatten

INTEGRATOR_TYPES = ('ERK', 'IRK', 'GNSF', 'DISCRETE', 'LIFTED_IRK')
COLLOCATION_TYPES = ('GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE', 'EXPLICIT_RUNGE_KUTTA')
COST_DISCRETIZATION_TYPES = ('EULER', 'INTEGRATOR')

class AcadosOcpOptions:
    """
    class containing the description of the solver options
    """
    def __init__(self):
        self.__hessian_approx = 'GAUSS_NEWTON'
        self.__integrator_type = 'ERK'
        self.__tf = None
        self.__N_horizon = None
        self.__nlp_solver_type = 'SQP'
        self.__nlp_solver_tol_stat = 1e-6
        self.__nlp_solver_tol_eq = 1e-6
        self.__nlp_solver_tol_ineq = 1e-6
        self.__nlp_solver_tol_comp = 1e-6
        self.__nlp_solver_max_iter = 100
        self.__nlp_solver_ext_qp_res = 0
        self.__nlp_solver_warm_start_first_qp = False
        self.__nlp_solver_warm_start_first_qp_from_nlp = False
        self.__nlp_solver_tol_min_step_norm = None
        self.__globalization = 'FIXED_STEP'
        self.__levenberg_marquardt = 0.0
        self.__collocation_type = 'GAUSS_LEGENDRE'
        self.__sim_method_num_stages = 4
        self.__sim_method_num_steps = 1
        self.__sim_method_newton_iter = 3
        self.__sim_method_newton_tol = 0.0
        self.__sim_method_jac_reuse = 0
        self.__shooting_nodes = None
        self.__time_steps = None
        self.__cost_scaling = None
        self.__Tsim = None
        self.__qp_solver = 'PARTIAL_CONDENSING_HPIPM'
        self.__qp_solver_tol_stat = None
        self.__qp_solver_tol_eq = None
        self.__qp_solver_tol_ineq = None
        self.__qp_solver_tol_comp = None
        self.__qp_solver_iter_max = 50
        self.__qp_solver_cond_N = None
        self.__qp_solver_cond_block_size = None
        self.__qp_solver_warm_start = 0
        self.__qp_solver_cond_ric_alg = 1
        self.__qp_solver_ric_alg = 1
        self.__qp_solver_mu0 = 0.0
        self.__qp_solver_t0_init = 2
        self.__tau_min = 0.0
        self.__solution_sens_qp_t_lam_min = 1e-9
        self.__rti_log_residuals = 0
        self.__rti_log_only_available_residuals = 0
        self.__print_level = 0
        self.__cost_discretization = 'EULER'
        self.__regularize_method = 'NO_REGULARIZE'
        self.__reg_epsilon = 1e-4
        self.__reg_max_cond_block = 1e7
        self.__reg_adaptive_eps = False
        self.__reg_min_epsilon = 1e-8
        self.__exact_hess_cost = 1
        self.__exact_hess_dyn = 1
        self.__exact_hess_constr = 1
        self.__eval_residual_at_max_iter = None
        self.__use_constraint_hessian_in_feas_qp = False
        self.__byrd_omojokon_slack_relaxation_factor = 1.00001
        self.__search_direction_mode = 'NOMINAL_QP'
        self.__allow_direction_mode_switch_to_nominal = True
        self.__fixed_hess = 0
        self.__globalization_funnel_init_increase_factor = 15.0
        self.__globalization_funnel_init_upper_bound = 1.0
        self.__globalization_funnel_sufficient_decrease_factor = 0.9
        self.__globalization_funnel_kappa = 0.9
        self.__globalization_funnel_fraction_switching_condition = 1e-3
        self.__globalization_funnel_initial_penalty_parameter = 1.0
        self.__globalization_funnel_use_merit_fun_only = False
        self.__globalization_fixed_step_length = 1.0
        self.__qpscaling_ub_max_abs_eig = 1e5
        self.__qpscaling_lb_norm_inf_grad_obj = 1e-4
        self.__qpscaling_scale_objective = "NO_OBJECTIVE_SCALING"
        self.__qpscaling_scale_constraints = "NO_CONSTRAINT_SCALING"

        self.__nlp_qp_tol_strategy = "FIXED_QP_TOL"
        self.__nlp_qp_tol_reduction_factor = 1e-1
        self.__nlp_qp_tol_safety_factor = 0.1
        self.__nlp_qp_tol_min_stat = 1e-9
        self.__nlp_qp_tol_min_eq = 1e-10
        self.__nlp_qp_tol_min_ineq = 1e-10
        self.__nlp_qp_tol_min_comp = 1e-11

        self.__ext_cost_num_hess = 0
        self.__globalization_use_SOC = 0
        self.__globalization_alpha_min = None
        self.__globalization_alpha_reduction = 0.7
        self.__globalization_line_search_use_sufficient_descent = 0
        self.__globalization_full_step_dual = None
        self.__globalization_eps_sufficient_descent = None
        self.__hpipm_mode = 'BALANCE'
        self.__with_solution_sens_wrt_params = False
        self.__with_value_sens_wrt_params = False
        self.__as_rti_iter = 1
        self.__as_rti_level = 4
        self.__with_adaptive_levenberg_marquardt = False
        self.__adaptive_levenberg_marquardt_lam = 5.0
        self.__adaptive_levenberg_marquardt_mu_min = 1e-16
        self.__adaptive_levenberg_marquardt_mu0 = 1e-3
        self.__adaptive_levenberg_marquardt_obj_scalar = 2.0
        self.__log_primal_step_norm: bool = False
        self.__log_dual_step_norm: bool = False
        self.__store_iterates: bool = False
        self.__timeout_max_time = 0.
        self.__timeout_heuristic = 'LAST'

        # TODO: move those out? they are more about generation than about the acados OCP solver.
        env = os.environ
        self.__ext_fun_compile_flags = '-O2' if 'ACADOS_EXT_FUN_COMPILE_FLAGS' not in env else env['ACADOS_EXT_FUN_COMPILE_FLAGS']
        self.__ext_fun_expand_constr = False
        self.__ext_fun_expand_cost = False
        self.__ext_fun_expand_precompute = False
        self.__ext_fun_expand_dyn = False
        self.__model_external_shared_lib_dir = None
        self.__model_external_shared_lib_name = None
        self.__custom_update_filename = ''
        self.__custom_update_header_filename = ''
        self.__custom_templates = []
        self.__custom_update_copy = True
        self.__num_threads_in_batch_solve: int = 1
        self.__with_batch_functionality: bool = False
        self.__with_anderson_acceleration: bool = False
        self.__anderson_activation_threshold: float = 1e1


    @property
    def qp_solver(self):
        """QP solver to be used in the NLP solver.
        String in ('PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP', 'PARTIAL_CONDENSING_CLARABEL', FULL_CONDENSING_DAQP').
        Default: 'PARTIAL_CONDENSING_HPIPM'.

        QP solver statuses are mapped to the acados status definitions.

        For HPIPM the status values are mapped as below
        | HPIPM status |   acados status   |
        |----------------------------------|
        |   SUCCESS   |  ACADOS_SUCCESS  0  |
        |    MAXIT    |  ACADOS_MAXITER  2  |
        |   MINSTEP   |  ACADOS_MINSTEP  3  |
        |     NAN     |    ACADOS_NAN    1  |
        |  INCONS_EQ  | ACADOS_INFEASIBLE 9 |
        |     ELSE    |  ACADOS_UNKNOWN -1  |

        For qpOASES the status values are mapped as below
        |       qpOASES status |        acados status         |
        |-----------------------------------------------------|
        |       SUCCESSFUL_RETURN       |  ACADOS_SUCCESS   0 |
        |      RET_MAX_NWSR_REACHED     |  ACADOS_MAXITER   2 |
        | RET_INIT_FAILED_UNBOUNDEDNESS | ACADOS_UNBOUNDED  6 |
        | RET_INIT_FAILED_INFEASIBILITY | ACADOS_INFEASIBLE 9 |
        |             ELSE              |  ACADOS_UNKNOWN  -1 |

        For DAQP the status values are mapped as below
        |     DAQP status     |       acados status      |
        |------------------------------------------------|
        |    EXIT_OPTIMAL     |     ACADOS_SUCCESS   0   |
        |  EXIT_SOFT_OPTIMAL  |     ACADOS_MAXITER   0   |
        | EXIT_EXIT_ITERLIMIT |     ACADOS_MAXITER   2   |
        |   EXIT_UNBOUNDED    |     ACADOS_UNBOUNDED 6   |
        |   EXIT_INFEASIBLE   |    ACADOS_INFEASIBLE 9   |
        |        ELSE         |     ACADOS_UNKNOWN  -1   |

        For QPDUNES the status values are mapped as below
        |             QPDUNES status            |       acados status         |
        |---------------------------------------------------------------------|
        |               QPDUNES_OK              |     ACADOS_SUCCESS    0     |
        | QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND   |     ACADOS_MAXITER    0     |
        | QPDUNES_ERR_ITERATION_LIMIT_REACHED   |     ACADOS_MAXITER    2     |
        |    QPDUNES_ERR_DIVISION_BY_ZERO       |    ACADOS_QP_FAILURE  4     |
        |   QPDUNES_ERR_STAGE_QP_INFEASIBLE     |    ACADOS_INFEASIBLE  9     |
        |                  ELSE                 |     ACADOS_UNKNOWN   -1     |

        """
        return self.__qp_solver

    @qp_solver.setter
    def qp_solver(self, qp_solver):
        qp_solvers = ('PARTIAL_CONDENSING_HPIPM', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM', \
                'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP', 'PARTIAL_CONDENSING_CLARABEL', \
                'FULL_CONDENSING_DAQP')
        if qp_solver in qp_solvers:
            self.__qp_solver = qp_solver
        else:
            raise ValueError('Invalid qp_solver value. Possible values are:\n\n' \
                    + ',\n'.join(qp_solvers) + '.\n\nYou have: ' + qp_solver + '.\n\n')


    @property
    def ext_fun_compile_flags(self):
        """
        String with compiler flags for external function compilation.
        Default: '-O2' if environment variable ACADOS_EXT_FUN_COMPILE_FLAGS is not set, else ACADOS_EXT_FUN_COMPILE_FLAGS is used as default.
        """
        return self.__ext_fun_compile_flags

    @ext_fun_compile_flags.setter
    def ext_fun_compile_flags(self, ext_fun_compile_flags):
        if isinstance(ext_fun_compile_flags, str):
            self.__ext_fun_compile_flags = ext_fun_compile_flags
        else:
            raise TypeError('Invalid ext_fun_compile_flags value, expected a string.\n')

    @property
    def ext_fun_expand_constr(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for constraint functions.
        Default: False
        """
        return self.__ext_fun_expand_constr

    @ext_fun_expand_constr.setter
    def ext_fun_expand_constr(self, ext_fun_expand_constr):
        if not isinstance(ext_fun_expand_constr, bool):
            raise TypeError('Invalid ext_fun_expand_constr value, expected bool.\n')
        self.__ext_fun_expand_constr = ext_fun_expand_constr

    @property
    def ext_fun_expand_cost(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for cost functions.
        Default: False
        """
        return self.__ext_fun_expand_cost

    @ext_fun_expand_cost.setter
    def ext_fun_expand_cost(self, ext_fun_expand_cost):
        if not isinstance(ext_fun_expand_cost, bool):
            raise TypeError('Invalid ext_fun_expand_cost value, expected bool.\n')
        self.__ext_fun_expand_cost = ext_fun_expand_cost

    @property
    def ext_fun_expand_dyn(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for dynamics functions.
        Default: False
        """
        return self.__ext_fun_expand_dyn

    @ext_fun_expand_dyn.setter
    def ext_fun_expand_dyn(self, ext_fun_expand_dyn):
        if not isinstance(ext_fun_expand_dyn, bool):
            raise TypeError('Invalid ext_fun_expand_dyn value, expected bool.\n')
        self.__ext_fun_expand_dyn = ext_fun_expand_dyn

    @property
    def ext_fun_expand_precompute(self):
        """
        Flag indicating whether CasADi.MX should be expanded to CasADi.SX before code generation for the precompute function.
        Default: False
        """
        return self.__ext_fun_expand_precompute

    @ext_fun_expand_precompute.setter
    def ext_fun_expand_precompute(self, ext_fun_expand_precompute):
        if not isinstance(ext_fun_expand_precompute, bool):
            raise TypeError('Invalid ext_fun_expand_precompute value, expected bool.\n')
        self.__ext_fun_expand_precompute = ext_fun_expand_precompute

    @property
    def custom_update_filename(self):
        """
        Filename of the custom C function to update solver data and parameters in between solver calls.
        Compare also `AcadosOcpOptions.custom_update_header_filename`.

        This file has to implement the functions
        int custom_update_init_function([model.name]_solver_capsule* capsule);
        int custom_update_function([model.name]_solver_capsule* capsule, double* data, int data_len);
        int custom_update_terminate_function([model.name]_solver_capsule* capsule);


        Default: ''.
        """
        return self.__custom_update_filename


    @custom_update_filename.setter
    def custom_update_filename(self, custom_update_filename):
        if isinstance(custom_update_filename, str):
            self.__custom_update_filename = custom_update_filename
        else:
            raise TypeError('Invalid custom_update_filename, expected a string.\n')

    @property
    def custom_templates(self):
        """
        List of tuples of the form:
        (input_filename, output_filename)

        Custom templates are render in OCP solver generation.

        Default: [].
        """
        return self.__custom_templates


    @custom_templates.setter
    def custom_templates(self, custom_templates):
        if not isinstance(custom_templates, list):
            raise TypeError('Invalid custom_templates, expected a list.\n')
        for tup in custom_templates:
            if not isinstance(tup, tuple):
                raise TypeError('Invalid custom_templates, should be list of tuples.\n')
            for s in tup:
                if not isinstance(s, str):
                    raise TypeError('Invalid custom_templates, should be list of tuples of strings.\n')
        self.__custom_templates = custom_templates

    @property
    def custom_update_header_filename(self):
        """
        Header filename of the custom C function to update solver data and parameters in between solver calls.

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

    @custom_update_header_filename.setter
    def custom_update_header_filename(self, custom_update_header_filename):
        if isinstance(custom_update_header_filename, str):
            self.__custom_update_header_filename = custom_update_header_filename
        else:
            raise TypeError('Invalid custom_update_header_filename, expected a string.\n')

    @property
    def custom_update_copy(self):
        """
        Boolean;
        If True, the custom update function files are copied into the `code_export_directory`.
        """
        return self.__custom_update_copy


    @custom_update_copy.setter
    def custom_update_copy(self, custom_update_copy):
        if isinstance(custom_update_copy, bool):
            self.__custom_update_copy = custom_update_copy
        else:
            raise TypeError('Invalid custom_update_copy, expected a bool.\n')

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

    @hpipm_mode.setter
    def hpipm_mode(self, hpipm_mode):
        hpipm_modes = ('BALANCE', 'SPEED_ABS', 'SPEED', 'ROBUST')
        if hpipm_mode in hpipm_modes:
            self.__hpipm_mode = hpipm_mode
        else:
            raise ValueError('Invalid hpipm_mode value. Possible values are:\n\n' \
                    + ',\n'.join(hpipm_modes) + '.\n\nYou have: ' + hpipm_mode + '.\n\n')

    @property
    def hessian_approx(self):
        """Hessian approximation.
        String in ('GAUSS_NEWTON', 'EXACT').
        Default: 'GAUSS_NEWTON'.
        """
        return self.__hessian_approx

    @hessian_approx.setter
    def hessian_approx(self, hessian_approx):
        hessian_approxs = ('GAUSS_NEWTON', 'EXACT')
        if hessian_approx in hessian_approxs:
            self.__hessian_approx = hessian_approx
        else:
            raise ValueError('Invalid hessian_approx value. Possible values are:\n\n' \
                    + ',\n'.join(hessian_approxs) + '.\n\nYou have: ' + hessian_approx + '.\n\n')

    @property
    def integrator_type(self):
        """
        Integrator type.
        String in ('ERK', 'IRK', 'GNSF', 'DISCRETE', 'LIFTED_IRK').
        Default: 'ERK'.
        """
        return self.__integrator_type

    @integrator_type.setter
    def integrator_type(self, integrator_type):
        if integrator_type in INTEGRATOR_TYPES:
            self.__integrator_type = integrator_type
        else:
            raise ValueError('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(INTEGRATOR_TYPES) + '.\n\nYou have: ' + integrator_type + '.\n\n')

    @property
    def nlp_solver_type(self):
        """NLP solver.
        String in ('SQP', 'SQP_RTI', 'DDP', 'SQP_WITH_FEASIBLE_QP').
        Default: 'SQP'.
        """
        return self.__nlp_solver_type

    @nlp_solver_type.setter
    def nlp_solver_type(self, nlp_solver_type):
        nlp_solver_types = ('SQP', 'SQP_RTI', 'DDP', 'SQP_WITH_FEASIBLE_QP')
        if nlp_solver_type in nlp_solver_types:
            self.__nlp_solver_type = nlp_solver_type
        else:
            raise ValueError('Invalid nlp_solver_type value. Possible values are:\n\n' \
                    + ',\n'.join(nlp_solver_types) + '.\n\nYou have: ' + nlp_solver_type + '.\n\n')

    @property
    def globalization(self):
        """Globalization type.
        String in ('FIXED_STEP', 'MERIT_BACKTRACKING', 'FUNNEL_L1PEN_LINESEARCH').
        Default: 'FIXED_STEP'.

        - FIXED_STEP: performs steps with a given step length, see option globalization_fixed_step_length
        - MERIT_BACKTRACKING: performs a merit function based backtracking line search following Following Leineweber1999, Section "3.5.1 Line Search Globalization"
        - FUNNEL_L1PEN_LINESEARCH: following "A Unified Funnel Restoration SQP Algorithm" by Kiessling et al.
            https://arxiv.org/pdf/2409.09208
            NOTE: preliminary implementation
        """
        return self.__globalization

    @globalization.setter
    def globalization(self, globalization):
        globalization_types = ('FUNNEL_L1PEN_LINESEARCH', 'MERIT_BACKTRACKING', 'FIXED_STEP')
        if globalization in globalization_types:
            self.__globalization = globalization
        else:
            raise ValueError('Invalid globalization value. Possible values are:\n\n' \
                    + ',\n'.join(globalization_types) + '.\n\nYou have: ' + globalization + '.\n\n')

    @property
    def collocation_type(self):
        """Collocation type: only relevant for implicit integrators
        -- string in {'GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE', 'EXPLICIT_RUNGE_KUTTA'}.

        Default: GAUSS_LEGENDRE.

        .. note:: GAUSS_LEGENDRE tableaus yield integration methods that are A-stable, but not L-stable and have order `2 * num_stages`,
        .. note:: GAUSS_RADAU_IIA tableaus yield integration methods that are L-stable and have order `2 * num_stages - 1`.
        .. note:: EXPLICIT_RUNGE_KUTTA tableaus can be used for comparisons of ERK and IRK to ensure correctness, but are only recommended with ERK for users.
        """
        return self.__collocation_type

    @collocation_type.setter
    def collocation_type(self, collocation_type):
        if collocation_type in COLLOCATION_TYPES:
            self.__collocation_type = collocation_type
        else:
            raise ValueError('Invalid collocation_type value. Possible values are:\n\n' \
                    + ',\n'.join(COLLOCATION_TYPES) + '.\n\nYou have: ' + collocation_type + '.\n\n')

    @property
    def regularize_method(self):
        """Regularization method for the Hessian.
        String in ('NO_REGULARIZE', 'MIRROR', 'PROJECT', 'PROJECT_REDUC_HESS', 'CONVEXIFY', 'GERSHGORIN_LEVENBERG_MARQUARDT').

        - MIRROR: performs eigenvalue decomposition H = V^T D V and sets D_ii = max(eps, abs(D_ii))
        - PROJECT: performs eigenvalue decomposition H = V^T D V and sets D_ii = max(eps, D_ii)
        - CONVEXIFY: Algorithm 6 from Verschueren2017, https://cdn.syscop.de/publications/Verschueren2017.pdf, experimental, might not be correct if inequality constraints are active.
        - PROJECT_REDUC_HESS: experimental, should make sure that the reduced Hessian is positive definite. Has to be used with qp_solver_ric_alg = 0 and qp_solver_cond_ric_alg = 0
        - GERSHGORIN_LEVENBERG_MARQUARDT: estimates the smallest eigenvalue of each Hessian block using Gershgorin circles and adds multiple of identity to each block, such that smallest eigenvalue after regularization is at least reg_epsilon

        Default: 'NO_REGULARIZE'.
        """
        return self.__regularize_method

    @regularize_method.setter
    def regularize_method(self, regularize_method):
        regularize_methods = ('NO_REGULARIZE', 'MIRROR', 'PROJECT', \
                                'PROJECT_REDUC_HESS', 'CONVEXIFY', 'GERSHGORIN_LEVENBERG_MARQUARDT')
        if regularize_method in regularize_methods:
            self.__regularize_method = regularize_method
        else:
            raise ValueError('Invalid regularize_method value. Possible values are:\n\n' \
                    + ',\n'.join(regularize_methods) + '.\n\nYou have: ' + regularize_method + '.\n\n')

    @property
    def globalization_fixed_step_length(self):
        """
        Fixed Newton step length, used if globalization == "FIXED_STEP"
        Type: float >= 0.
        Default: 1.0.
        """
        return self.__globalization_fixed_step_length

    @globalization_fixed_step_length.setter
    def globalization_fixed_step_length(self, globalization_fixed_step_length):
        if isinstance(globalization_fixed_step_length, float) and globalization_fixed_step_length >= 0 or globalization_fixed_step_length <= 1.0:
            self.__globalization_fixed_step_length = globalization_fixed_step_length
        else:
            raise ValueError('Invalid globalization_fixed_step_length value. globalization_fixed_step_length must be a float in [0, 1].')

    @property
    def qpscaling_ub_max_abs_eig(self):
        """
        Maximum upper bound for eigenvalues of Hessian after QP scaling.
        Type: float >= 0.
        Default: 1e5.
        """
        return self.__qpscaling_ub_max_abs_eig

    @qpscaling_ub_max_abs_eig.setter
    def qpscaling_ub_max_abs_eig(self, qpscaling_ub_max_abs_eig):
        if isinstance(qpscaling_ub_max_abs_eig, float) and qpscaling_ub_max_abs_eig >= 0.:
            self.__qpscaling_ub_max_abs_eig = qpscaling_ub_max_abs_eig
        else:
            raise ValueError('Invalid qpscaling_ub_max_abs_eig value. qpscaling_ub_max_abs_eig must be a positive float.')

    @property
    def qpscaling_lb_norm_inf_grad_obj(self):
        """
        Minimum allowed lower bound for inf norm of qp gradient in QP scaling.
        Is attempted to be respected, respecting qpscaling_ub_max_abs_eig is prioritized.
        Type: float >= 0.
        Default: 1e-4.
        """
        return self.__qpscaling_lb_norm_inf_grad_obj

    @qpscaling_lb_norm_inf_grad_obj.setter
    def qpscaling_lb_norm_inf_grad_obj(self, qpscaling_lb_norm_inf_grad_obj):
        if isinstance(qpscaling_lb_norm_inf_grad_obj, float) and qpscaling_lb_norm_inf_grad_obj >= 0.:
            self.__qpscaling_lb_norm_inf_grad_obj = qpscaling_lb_norm_inf_grad_obj
        else:
            raise ValueError('Invalid qpscaling_lb_norm_inf_grad_obj value. qpscaling_lb_norm_inf_grad_obj must be a positive float.')

    @property
    def qpscaling_scale_objective(self):
        """
        String in ["NO_OBJECTIVE_SCALING", "OBJECTIVE_GERSHGORIN"]
        Default: "NO_OBJECTIVE_SCALING".

        - NO_OBJECTIVE_SCALING: no scaling of the objective
        - OBJECTIVE_GERSHGORIN: estimate max. abs. eigenvalue using Gershgorin circles as `max_abs_eig`, then sets the objective scaling factor as `obj_factor = min(1.0, qpscaling_ub_max_abs_eig/max_abs_eig)`
        """
        return self.__qpscaling_scale_objective

    @qpscaling_scale_objective.setter
    def qpscaling_scale_objective(self, qpscaling_scale_objective):
        qpscaling_scale_objective_types = ["NO_OBJECTIVE_SCALING", "OBJECTIVE_GERSHGORIN"]
        if not qpscaling_scale_objective in qpscaling_scale_objective_types:
            raise ValueError(f'Invalid qpscaling_scale_objective value. Must be in {qpscaling_scale_objective_types}, got {qpscaling_scale_objective}.')
        self.__qpscaling_scale_objective = qpscaling_scale_objective

    @property
    def qpscaling_scale_constraints(self):
        """
        String in ["NO_CONSTRAINT_SCALING", "INF_NORM"]
        Default: "NO_CONSTRAINT_SCALING".

        - NO_CONSTRAINT_SCALING: no scaling of the constraints
        - INF_NORM: scales each constraint except simple bounds by factor `1.0 / max(inf_norm_coeff, inf_norm_constraint_bound)`, such that the inf-norm of the constraint coefficients and bounds is <= 1.0.
        Slack penalties are adjusted accordingly to get an equivalent solution.
        First, the cost is scaled, then the constraints.
        """
        return self.__qpscaling_scale_constraints

    @qpscaling_scale_constraints.setter
    def qpscaling_scale_constraints(self, qpscaling_scale_constraints):
        qpscaling_scale_constraints_types = ["NO_CONSTRAINT_SCALING", "INF_NORM"]
        if not qpscaling_scale_constraints in qpscaling_scale_constraints_types:
            raise ValueError(f'Invalid qpscaling_scale_constraints value. Must be in {qpscaling_scale_constraints_types}, got {qpscaling_scale_constraints}.')
        self.__qpscaling_scale_constraints = qpscaling_scale_constraints

    @property
    def nlp_qp_tol_strategy(self):
        """
        Strategy for setting the QP tolerances in the NLP solver.
        String in ["ADAPTIVE_CURRENT_RES_JOINT", "ADAPTIVE_QPSCALING", "FIXED_QP_TOL"]

        - FIXED_QP_TOL: uses the fixed QP solver tolerances set by the properties `qp_solver_tol_stat`, `qp_solver_tol_eq`, `qp_solver_tol_ineq`, `qp_solver_tol_comp`, only this was implemented in acados <= v0.5.0.

        - ADAPTIVE_CURRENT_RES_JOINT: uses the current NLP residuals to set the QP tolerances in a joint manner.
        The QP tolerances are set as follows:
            1) `tmp_tol_* = MIN(nlp_qp_tol_reduction_factor * inf_norm_res_*, 1e-2)`
            2) `joint_tol = MAX(tmp_tol_* for all * in ['stat', 'eq', 'ineq', 'comp'])`
            3) `tol_* = MAX(joint_tol, nlp_qp_tol_safety_factor * nlp_solver_tol_*)`

        - ADAPTIVE_QPSCALING: adapts the QP tolerances based on the QP scaling factors, to make NLP residuals converge to desired tolerances, if it can be achieved.
        The QP tolerances are set as follows:
            1) `qp_tol_stat = nlp_qp_tol_safety_factor * nlp_solver_tol_stat * MIN(objective_scaling_factor, min_constraint_scaling);`
            2) `qp_tol_eq = nlp_qp_tol_safety_factor * nlp_solver_tol_eq`
            3) `qp_tol_ineq = nlp_qp_tol_safety_factor * nlp_solver_tol_ineq * min_constraint_scaling`
            4) `qp_tol_comp = nlp_qp_tol_safety_factor * nlp_solver_tol_comp * min_constraint_scaling`
            5) cap all QP tolerances to a minimum of `nlp_qp_tol_min_*`.

        Default: "FIXED_QP_TOL".
        """
        return self.__nlp_qp_tol_strategy

    @nlp_qp_tol_strategy.setter
    def nlp_qp_tol_strategy(self, nlp_qp_tol_strategy):
        nlp_qp_tol_strategy_types = ["ADAPTIVE_CURRENT_RES_JOINT", "ADAPTIVE_QPSCALING", "FIXED_QP_TOL"]
        if not nlp_qp_tol_strategy in nlp_qp_tol_strategy_types:
            raise ValueError(f'Invalid nlp_qp_tol_strategy value. Must be in {nlp_qp_tol_strategy_types}, got {nlp_qp_tol_strategy}.')
        self.__nlp_qp_tol_strategy = nlp_qp_tol_strategy

    @property
    def nlp_qp_tol_reduction_factor(self):
        """
        Factor by which the QP tolerance is smaller compared to the NLP residuals when using the ADAPTIVE_CURRENT_RES_JOINT strategy.
        Default: 1e-1.
        """
        return self.__nlp_qp_tol_reduction_factor

    @nlp_qp_tol_reduction_factor.setter
    def nlp_qp_tol_reduction_factor(self, nlp_qp_tol_reduction_factor):
        if not isinstance(nlp_qp_tol_reduction_factor, float) or nlp_qp_tol_reduction_factor < 0.0 or nlp_qp_tol_reduction_factor > 1.0:
            raise ValueError(f'Invalid nlp_qp_tol_reduction_factor value. Must be in [0, 1], got {nlp_qp_tol_reduction_factor}.')
        self.__nlp_qp_tol_reduction_factor = nlp_qp_tol_reduction_factor

    @property
    def nlp_qp_tol_safety_factor(self):
        """
        Safety factor for the QP tolerances.
        Used to ensure qp_tol* = nlp_qp_tol_safety_factor * nlp_solver_tol_* when approaching the NLP solution.
        Often QPs should be solved to a higher accuracy than the NLP solver tolerances, to ensure convergence of the NLP solver.
        Used in the ADAPTIVE_CURRENT_RES_JOINT, ADAPTIVE_QPSCALING strategies.
        Type: float in [0, 1].
        Default: 0.1.
        """
        return self.__nlp_qp_tol_safety_factor

    @nlp_qp_tol_safety_factor.setter
    def nlp_qp_tol_safety_factor(self, nlp_qp_tol_safety_factor):
        if not isinstance(nlp_qp_tol_safety_factor, float) or nlp_qp_tol_safety_factor < 0.0 or nlp_qp_tol_safety_factor > 1.0:
            raise ValueError(f'Invalid nlp_qp_tol_safety_factor value. Must be in [0, 1], got {nlp_qp_tol_safety_factor}.')
        self.__nlp_qp_tol_safety_factor = nlp_qp_tol_safety_factor

    @property
    def nlp_qp_tol_min_stat(self):
        """
        Minimum value to be set in the QP solver stationarity tolerance by `nlp_qp_tol_strategy`, used in `ADAPTIVE_QPSCALING`.
        Type: float > 0.
        Default: 1e-9.
        """
        return self.__nlp_qp_tol_min_stat

    @property
    def nlp_qp_tol_min_eq(self):
        """
        Minimum value to be set in the QP solver equality tolerance by `nlp_qp_tol_strategy`, used in `ADAPTIVE_QPSCALING`.
        Type: float > 0.
        Default: 1e-10.
        """
        return self.__nlp_qp_tol_min_eq

    @property
    def nlp_qp_tol_min_ineq(self):
        """
        Minimum value to be set in the QP solver inequality tolerance by `nlp_qp_tol_strategy`, used in `ADAPTIVE_QPSCALING`.
        Type: float > 0.
        Default: 1e-10.
        """
        return self.__nlp_qp_tol_min_ineq

    @property
    def nlp_qp_tol_min_comp(self):
        """
        Minimum value to be set in the QP solver complementarity tolerance by `nlp_qp_tol_strategy`, used in `ADAPTIVE_QPSCALING`.
        Type: float > 0.
        Default: 1e-11.
        """
        return self.__nlp_qp_tol_min_comp

    @property
    @deprecated(version="0.4.0", reason="Use globalization_fixed_step_length instead.")
    def nlp_solver_step_length(self):
        """
        This option is deprecated and has new name: globalization_fixed_step_length
        """
        return self.__globalization_fixed_step_length

    @nlp_solver_step_length.setter
    @deprecated(version="0.4.0", reason="Use globalization_fixed_step_length instead.")
    def nlp_solver_step_length(self, nlp_solver_step_length):
        self.globalization_fixed_step_length = nlp_solver_step_length

    @property
    def nlp_solver_warm_start_first_qp(self):
        """
        Flag indicating whether the first QP in an NLP solve should be warm started.
        If the warm start should be done using the NLP iterate, see property nlp_solver_warm_start_first_qp_from_nlp.

        The warm start level of the QP solver is controlled by the property qp_solver_warm_start.

        Type: bool.
        Default: False.
        """
        return self.__nlp_solver_warm_start_first_qp

    @nlp_solver_warm_start_first_qp.setter
    def nlp_solver_warm_start_first_qp(self, nlp_solver_warm_start_first_qp):
        if isinstance(nlp_solver_warm_start_first_qp, bool):
            self.__nlp_solver_warm_start_first_qp = nlp_solver_warm_start_first_qp
        else:
            raise TypeError('Invalid nlp_solver_warm_start_first_qp value. Expected bool.')

    @property
    def nlp_solver_warm_start_first_qp_from_nlp(self):
        """
        Only relevant if `nlp_solver_warm_start_first_qp` is True.
        If True first QP will be initialized using values from NLP iterate, otherwise from previous QP solution.

        Note: for now only works with HPIPM and partial condensing with N = qp_solver_partial_cond_N
        Type: bool.
        Default: False.
        """
        return self.__nlp_solver_warm_start_first_qp_from_nlp


    @nlp_solver_warm_start_first_qp_from_nlp.setter
    def nlp_solver_warm_start_first_qp_from_nlp(self, nlp_solver_warm_start_first_qp_from_nlp):
        if not isinstance(nlp_solver_warm_start_first_qp_from_nlp, bool):
            raise TypeError('Invalid nlp_solver_warm_start_first_qp_from_nlp value. Expected bool.')
        self.__nlp_solver_warm_start_first_qp_from_nlp = nlp_solver_warm_start_first_qp_from_nlp

    @property
    def levenberg_marquardt(self):
        """
        Factor for LM regularization.
        Type: float >= 0
        Default: 0.0.
        """
        return self.__levenberg_marquardt

    @levenberg_marquardt.setter
    def levenberg_marquardt(self, levenberg_marquardt):
        if isinstance(levenberg_marquardt, float) and levenberg_marquardt >= 0:
            self.__levenberg_marquardt = levenberg_marquardt
        else:
            raise ValueError('Invalid levenberg_marquardt value. levenberg_marquardt must be a positive float.')

    @property
    def sim_method_num_stages(self):
        """
        Number of stages in the integrator.
        Type: int > 0 or ndarray of ints > 0 of shape (N,).
        Default: 4
        """
        return self.__sim_method_num_stages

    @sim_method_num_stages.setter
    def sim_method_num_stages(self, sim_method_num_stages):
        # NOTE: checks in make_consistent
        self.__sim_method_num_stages = sim_method_num_stages

    @property
    def sim_method_num_steps(self):
        """
        Number of steps in the integrator.
        Type: int > 0 or ndarray of ints > 0 of shape (N,).
        Default: 1
        """
        return self.__sim_method_num_steps

    @sim_method_num_steps.setter
    def sim_method_num_steps(self, sim_method_num_steps):
        # NOTE: checks in make_consistent
        self.__sim_method_num_steps = sim_method_num_steps


    @property
    def sim_method_newton_iter(self):
        """
        Number of Newton iterations in implicit integrators.
        Type: int > 0
        Default: 3
        """
        return self.__sim_method_newton_iter

    @sim_method_newton_iter.setter
    def sim_method_newton_iter(self, sim_method_newton_iter):

        if isinstance(sim_method_newton_iter, int):
            self.__sim_method_newton_iter = sim_method_newton_iter
        else:
            raise ValueError('Invalid sim_method_newton_iter value. sim_method_newton_iter must be an integer.')

    @property
    def sim_method_newton_tol(self):
        """
        Tolerance of Newton system in implicit integrators.
        This option is not implemented for LIFTED_IRK
        Type: float: 0.0 means not used
        Default: 0.0
        """
        return self.__sim_method_newton_tol

    @sim_method_newton_tol.setter
    def sim_method_newton_tol(self, sim_method_newton_tol):
        if isinstance(sim_method_newton_tol, float) and sim_method_newton_tol > 0:
            self.__sim_method_newton_tol = sim_method_newton_tol
        else:
            raise ValueError('Invalid sim_method_newton_tol value. sim_method_newton_tol must be a positive float.')

    @property
    def sim_method_jac_reuse(self):
        """
        Integer determining if jacobians are reused within integrator or ndarray of ints > 0 of shape (N,).
        0: False (no reuse); 1: True (reuse)
        Default: 0
        """
        return self.__sim_method_jac_reuse

    @sim_method_jac_reuse.setter
    def sim_method_jac_reuse(self, sim_method_jac_reuse):
        self.__sim_method_jac_reuse = sim_method_jac_reuse

    @property
    def qp_solver_tol_stat(self):
        """
        QP solver stationarity tolerance.
        Used if nlp_qp_tol_strategy == "FIXED_QP_TOL".
        Default: :code:`None`
        """
        return self.__qp_solver_tol_stat

    @qp_solver_tol_stat.setter
    def qp_solver_tol_stat(self, qp_solver_tol_stat):
        if isinstance(qp_solver_tol_stat, float) and qp_solver_tol_stat > 0:
            self.__qp_solver_tol_stat = qp_solver_tol_stat
        else:
            raise ValueError('Invalid qp_solver_tol_stat value. qp_solver_tol_stat must be a positive float.')

    @property
    def qp_solver_tol_eq(self):
        """
        QP solver equality tolerance.
        Used if nlp_qp_tol_strategy == "FIXED_QP_TOL".
        Default: :code:`None`
        """
        return self.__qp_solver_tol_eq

    @qp_solver_tol_eq.setter
    def qp_solver_tol_eq(self, qp_solver_tol_eq):
        if isinstance(qp_solver_tol_eq, float) and qp_solver_tol_eq > 0:
            self.__qp_solver_tol_eq = qp_solver_tol_eq
        else:
            raise ValueError('Invalid qp_solver_tol_eq value. qp_solver_tol_eq must be a positive float.')

    @property
    def qp_solver_tol_ineq(self):
        """
        QP solver inequality.
        Used if nlp_qp_tol_strategy == "FIXED_QP_TOL".
        Default: :code:`None`
        """
        return self.__qp_solver_tol_ineq

    @qp_solver_tol_ineq.setter
    def qp_solver_tol_ineq(self, qp_solver_tol_ineq):
        if isinstance(qp_solver_tol_ineq, float) and qp_solver_tol_ineq > 0:
            self.__qp_solver_tol_ineq = qp_solver_tol_ineq
        else:
            raise ValueError('Invalid qp_solver_tol_ineq value. qp_solver_tol_ineq must be a positive float.')

    @property
    def qp_solver_tol_comp(self):
        """
        QP solver complementarity.
        Used if nlp_qp_tol_strategy == "FIXED_QP_TOL".
        Default: :code:`None`
        """
        return self.__qp_solver_tol_comp

    @qp_solver_tol_comp.setter
    def qp_solver_tol_comp(self, qp_solver_tol_comp):
        if isinstance(qp_solver_tol_comp, float) and qp_solver_tol_comp > 0:
            self.__qp_solver_tol_comp = qp_solver_tol_comp
        else:
            raise ValueError('Invalid qp_solver_tol_comp value. qp_solver_tol_comp must be a positive float.')

    @property
    def qp_solver_cond_N(self):
        """QP solver: New horizon after partial condensing.
        Set to N by default -> no condensing."""
        return self.__qp_solver_cond_N

    @qp_solver_cond_N.setter
    def qp_solver_cond_N(self, qp_solver_cond_N):
        if isinstance(qp_solver_cond_N, int) and qp_solver_cond_N >= 0:
            self.__qp_solver_cond_N = qp_solver_cond_N
        else:
            raise ValueError('Invalid qp_solver_cond_N value. qp_solver_cond_N must be a positive int.')

    @property
    def qp_solver_cond_block_size(self):
        """QP solver: list of integers of length qp_solver_cond_N + 1
        Denotes how many blocks of the original OCP are lumped together into one in partial condensing.
        Note that the last entry is the number of blocks that are condensed into the terminal cost of the partially condensed QP.
        Default: None -> compute even block size distribution based on qp_solver_cond_N
        """
        return self.__qp_solver_cond_block_size


    @qp_solver_cond_block_size.setter
    def qp_solver_cond_block_size(self, qp_solver_cond_block_size):
        if not isinstance(qp_solver_cond_block_size, list):
            raise ValueError('Invalid qp_solver_cond_block_size value. qp_solver_cond_block_size must be a list of nonnegative integers.')
        for i in qp_solver_cond_block_size:
            if not isinstance(i, int) or not i >= 0:
                raise ValueError('Invalid qp_solver_cond_block_size value. qp_solver_cond_block_size must be a list of nonnegative integers.')
        self.__qp_solver_cond_block_size = qp_solver_cond_block_size

    @property
    def qp_solver_warm_start(self):
        """
        Controls the QP solver warm start level.
        The very first QP in an NLP solve is by default not warm started, i.e. the QP warm start level for the first QP solve is set to 0.
        To warm start also the first QP, set nlp_solver_warm_start_first_qp.
        Also see nlp_solver_warm_start_first_qp_from_nlp.

        What warm/hot start means in detail is dependend on the QP solver being used.
        0: no warm start; 1: warm start; 2: hot start.
        Default: 0
        """
        return self.__qp_solver_warm_start

    @qp_solver_warm_start.setter
    def qp_solver_warm_start(self, qp_solver_warm_start):
        if qp_solver_warm_start in [0, 1, 2, 3]:
            self.__qp_solver_warm_start = qp_solver_warm_start
        else:
            raise ValueError('Invalid qp_solver_warm_start value. qp_solver_warm_start must be 0 or 1 or 2 or 3.')

    @property
    def qp_solver_cond_ric_alg(self):
        """
        QP solver: Determines which algorithm is used in HPIPM condensing.
        0: dont factorize hessian in the condensing; 1: factorize.
        Default: 1
        """
        return self.__qp_solver_cond_ric_alg

    @qp_solver_cond_ric_alg.setter
    def qp_solver_cond_ric_alg(self, qp_solver_cond_ric_alg):
        if qp_solver_cond_ric_alg in [0, 1]:
            self.__qp_solver_cond_ric_alg = qp_solver_cond_ric_alg
        else:
            raise ValueError(f'Invalid qp_solver_cond_ric_alg value. qp_solver_cond_ric_alg must be in [0, 1], got {qp_solver_cond_ric_alg}.')


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

    @qp_solver_ric_alg.setter
    def qp_solver_ric_alg(self, qp_solver_ric_alg):
        if qp_solver_ric_alg in [0, 1]:
            self.__qp_solver_ric_alg = qp_solver_ric_alg
        else:
            raise ValueError(f'Invalid qp_solver_ric_alg value. qp_solver_ric_alg must be in [0, 1], got {qp_solver_ric_alg}.')

    @property
    def qp_solver_mu0(self):
        """
        For HPIPM QP solver: Initial value for the barrier parameter.
        If 0, the default value according to hpipm_mode is used.

        Default: 0
        """
        return self.__qp_solver_mu0

    @qp_solver_mu0.setter
    def qp_solver_mu0(self, qp_solver_mu0):
        if isinstance(qp_solver_mu0, float) and qp_solver_mu0 >= 0:
            self.__qp_solver_mu0 = qp_solver_mu0
        else:
            raise ValueError('Invalid qp_solver_mu0 value. qp_solver_mu0 must be a positive float.')

    @property
    def qp_solver_t0_init(self):
        """
        For HPIPM QP solver: Initialization scheme of lambda and t slacks within HPIPM.
        0: initialize with sqrt(mu0)
        1: initialize with 1.0
        2: heuristic for primal feasibility

        When using larger value for tau_min, it is beneficial to not use 2, as the initialization of (t, lambda) might be too far off from the central path and prevent convergence.

        Type: int > 0
        Default: 2
        """
        return self.__qp_solver_t0_init

    @qp_solver_t0_init.setter
    def qp_solver_t0_init(self, qp_solver_t0_init):
        if qp_solver_t0_init in [0, 1, 2]:
            self.__qp_solver_t0_init = qp_solver_t0_init
        else:
            raise ValueError('Invalid qp_solver_t0_init value. Must be in [0, 1, 2].')

    @property
    def tau_min(self):
        """
        Minimum value for the barrier parameter tau.
        Relevant if an interior point method is used as a (sub)solver, right now this is only HPIPM.
        If no interior point method is used, this is set to 0.

        Default: 0
        """
        return self.__tau_min

    @tau_min.setter
    def tau_min(self, tau_min):
        if isinstance(tau_min, float) and tau_min >= 0:
            self.__tau_min = tau_min
        else:
            raise ValueError('Invalid tau_min value. tau_min must be a positive float.')

    @property
    def solution_sens_qp_t_lam_min(self):
        """
        When computing the solution sensitivities using the function `setup_qp_matrices_and_factorize()`, this value is used to clip the values lambda and t slack values of the QP iterate before factorization.

        Default: 1e-9
        """
        return self.__solution_sens_qp_t_lam_min

    @solution_sens_qp_t_lam_min.setter
    def solution_sens_qp_t_lam_min(self, solution_sens_qp_t_lam_min):
        if isinstance(solution_sens_qp_t_lam_min, float) and solution_sens_qp_t_lam_min >= 0:
            self.__solution_sens_qp_t_lam_min = solution_sens_qp_t_lam_min
        else:
            raise ValueError('Invalid solution_sens_qp_t_lam_min value. solution_sens_qp_t_lam_min must be a nonnegative float.')

    @property
    def qp_solver_iter_max(self):
        """
        QP solver: maximum number of iterations.
        Type: int > 0
        Default: 50
        """
        return self.__qp_solver_iter_max


    @qp_solver_iter_max.setter
    def qp_solver_iter_max(self, qp_solver_iter_max):
        if isinstance(qp_solver_iter_max, int) and qp_solver_iter_max >= 0:
            self.__qp_solver_iter_max = qp_solver_iter_max
        else:
            raise ValueError('Invalid qp_solver_iter_max value. qp_solver_iter_max must be a positive int.')

    @property
    def as_rti_iter(self):
        """
        Maximum number of iterations in the advanced-step real-time iteration.
        Default: 1
        """
        return self.__as_rti_iter

    @as_rti_iter.setter
    def as_rti_iter(self, as_rti_iter):
        if isinstance(as_rti_iter, int) and as_rti_iter >= 0:
            self.__as_rti_iter = as_rti_iter
        else:
            raise ValueError('Invalid as_rti_iter value. as_rti_iter must be a nonnegative int.')

    @property
    def with_anderson_acceleration(self):
        """
        Determines if algorithm uses Anderson accelerations.
        Only depth one is supported.
        Anderson accelerations are performed whenever the infinity norm of the KKT residual is < `anderson_activation_threshold `.
        Only supported for globalization == 'FIXED_STEP'.

        Type: bool
        Default: False
        """
        return self.__with_anderson_acceleration

    @with_anderson_acceleration.setter
    def with_anderson_acceleration(self, with_anderson_acceleration):
        if not isinstance(with_anderson_acceleration, bool):
            raise TypeError('Invalid with_anderson_acceleration value, must be bool.')
        self.__with_anderson_acceleration = with_anderson_acceleration

    @property
    def anderson_activation_threshold(self):
        """
        Only relevant if with_anderson_acceleration == True.
        Anderson accelerations are performed whenever the infinity norm of the KKT residual is < `anderson_activation_threshold `.

        If the KKT residual norm is larger than `anderson_activation_threshold `, no Anderson acceleration is performed.

        In the language of [Pollock2021, Sec. 5.1]*, this corresponds to specifying an "initial regime", consisting of iterates with KKT residual norm > `anderson_activation_threshold ` and an (pre-)asymptotic regime, where the residual norm is <= `anderson_activation_threshold`.
        In the initial regime, no Anderson acceleration is performed, i.e. depth $m=0$.
        In the (pre-)asymptotic regime, Anderson acceleration with depth $m=1$ is performed.

        *[Pollock2021] Anderson acceleration for contractive and noncontractive operators, Sara Pollock, IMA Journal of Numerical Analysis, 2021

        Type: float
        Default: 1e1
        """
        return self.__anderson_activation_threshold


    @anderson_activation_threshold.setter
    def anderson_activation_threshold(self, anderson_activation_threshold):
        if not isinstance(anderson_activation_threshold, float):
            raise TypeError('Invalid anderson_activation_threshold value, must be float.')
        self.__anderson_activation_threshold = anderson_activation_threshold

    @property
    def as_rti_level(self):
        """
        Level of the advanced-step real-time iteration.

        LEVEL-A: 0
        LEVEL-B: 1
        LEVEL-C: 2
        LEVEL-D: 3
        STANDARD_RTI: 4

        Default: 4
        """
        return self.__as_rti_level

    @as_rti_level.setter
    def as_rti_level(self, as_rti_level):
        if as_rti_level in [0, 1, 2, 3, 4]:
            self.__as_rti_level = as_rti_level
        else:
            raise ValueError('Invalid as_rti_level value must be in [0, 1, 2, 3, 4].')


    @property
    def with_adaptive_levenberg_marquardt(self):
        """
        So far: Only relevant for DDP
        Flag indicating whether adaptive levenberg marquardt should be used.
        This is useful for NLS or LS problem where the residual goes to zero since
        quadratic local convergence can be achieved.
        type: bool
        """
        return self.__with_adaptive_levenberg_marquardt

    @with_adaptive_levenberg_marquardt.setter
    def with_adaptive_levenberg_marquardt(self, with_adaptive_levenberg_marquardt):
        if isinstance(with_adaptive_levenberg_marquardt, bool):
            self.__with_adaptive_levenberg_marquardt = with_adaptive_levenberg_marquardt
        else:
            raise TypeError('Invalid with_adaptive_levenberg_marquardt value. Expected bool.')

    @property
    def adaptive_levenberg_marquardt_lam(self):
        """
        So far: Only relevant for DDP
        Flag defining the value of lambda in the adaptive levenberg_marquardt.
        Must be > 1.
        type: float
        """
        return self.__adaptive_levenberg_marquardt_lam

    @adaptive_levenberg_marquardt_lam.setter
    def adaptive_levenberg_marquardt_lam(self, adaptive_levenberg_marquardt_lam):
        if isinstance(adaptive_levenberg_marquardt_lam, float) and adaptive_levenberg_marquardt_lam >= 1.0:
            self.__adaptive_levenberg_marquardt_lam = adaptive_levenberg_marquardt_lam
        else:
            raise ValueError('Invalid adaptive_levenberg_marquardt_lam value. adaptive_levenberg_marquardt_lam must be a float greater 1.0.')

    @property
    def adaptive_levenberg_marquardt_mu_min(self):
        """
        So far: Only relevant for DDP
        Flag defining the value of mu_min in the adaptive levenberg_marquardt.
        Must be > 0.
        type: float
        """
        return self.__adaptive_levenberg_marquardt_mu_min

    @adaptive_levenberg_marquardt_mu_min.setter
    def adaptive_levenberg_marquardt_mu_min(self, adaptive_levenberg_marquardt_mu_min):
        if isinstance(adaptive_levenberg_marquardt_mu_min, float) and adaptive_levenberg_marquardt_mu_min >= 0.0:
            self.__adaptive_levenberg_marquardt_mu_min = adaptive_levenberg_marquardt_mu_min
        else:
            raise ValueError('Invalid adaptive_levenberg_marquardt_mu_min value. adaptive_levenberg_marquardt_mu_min must be a positive float.')

    @property
    def adaptive_levenberg_marquardt_mu0(self):
        """
        So far: Only relevant for DDP
        Flag defining the value of mu0 in the adaptive levenberg_marquardt.
        Must be > 0.
        type: float
        """
        return self.__adaptive_levenberg_marquardt_mu0

    @adaptive_levenberg_marquardt_mu0.setter
    def adaptive_levenberg_marquardt_mu0(self, adaptive_levenberg_marquardt_mu0):
        if isinstance(adaptive_levenberg_marquardt_mu0, float) and adaptive_levenberg_marquardt_mu0 >= 0.0:
            self.__adaptive_levenberg_marquardt_mu0 = adaptive_levenberg_marquardt_mu0
        else:
            raise ValueError('Invalid adaptive_levenberg_marquardt_mu0 value. adaptive_levenberg_marquardt_mu0 must be a positive float.')

    @property
    def adaptive_levenberg_marquardt_obj_scalar(self):
        """
        So far: Only relevant for DDP
        Flag defining the value of the scalar that is multiplied with the NLP
        objective function in the adaptive levenberg_marquardt.
        Must be > 0.
        Default: 2.0
        type: float
        """
        return self.__adaptive_levenberg_marquardt_obj_scalar

    @adaptive_levenberg_marquardt_obj_scalar.setter
    def adaptive_levenberg_marquardt_obj_scalar(self, adaptive_levenberg_marquardt_obj_scalar):
        if isinstance(adaptive_levenberg_marquardt_obj_scalar, float) and adaptive_levenberg_marquardt_obj_scalar >= 0.0:
            self.__adaptive_levenberg_marquardt_obj_scalar = adaptive_levenberg_marquardt_obj_scalar
        else:
            raise ValueError('Invalid adaptive_levenberg_marquardt_obj_scalar value. adaptive_levenberg_marquardt_obj_scalar must be a positive float.')

    @property
    def log_primal_step_norm(self):
        """
        Flag indicating whether the max norm of the primal steps should be logged.
        This is implemented only for solver types `SQP`, `SQP_WITH_FEASIBLE_QP`.
        Default: False
        """
        return self.__log_primal_step_norm

    @log_primal_step_norm.setter
    def log_primal_step_norm(self, val):
        if not isinstance(val, bool):
            raise TypeError('Invalid log_primal_step_norm value. Expected bool.')
        self.__log_primal_step_norm = val

    @property
    def log_dual_step_norm(self):
        """
        Flag indicating whether the max norm of the dual steps should be logged.
        This is implemented only for solver types `SQP`, `SQP_WITH_FEASIBLE_QP`.
        Default: False
        """
        return self.__log_dual_step_norm

    @log_dual_step_norm.setter
    def log_dual_step_norm(self, val):
        if not isinstance(val, bool):
            raise TypeError('Invalid log_dual_step_norm value. Expected bool.')
        self.__log_dual_step_norm = val

    @property
    def store_iterates(self):
        """
        Flag indicating whether the intermediate primal-dual iterates should be stored.
        This is implemented only for solver types `SQP` and `DDP`.
        Default: False
        """
        return self.__store_iterates

    @store_iterates.setter
    def store_iterates(self, val):
        if isinstance(val, bool):
            self.__store_iterates = val
        else:
            raise TypeError('Invalid store_iterates value. Expected bool.')

    @property
    def timeout_max_time(self):
        """
        Maximum time before solver timeout. If 0, there is no timeout.
        A timeout is triggered if the condition
        `current_time_tot + predicted_per_iteration_time > timeout_max_time`
        is satisfied at the end of an SQP iteration.
        The value of `predicted_per_iteration_time` is estimated using `timeout_heuristic`.
        Currently implemented for SQP only.
        Default: 0.
        """
        return self.__timeout_max_time

    @timeout_max_time.setter
    def timeout_max_time(self, val):
        if isinstance(val, float) and val >= 0:
            self.__timeout_max_time = val
        else:
            raise ValueError('Invalid timeout_max_time value. Expected nonnegative float.')

    @property
    def timeout_heuristic(self):
        """
        Heuristic to be used for predicting the runtime of the next SQP iteration, cf. `timeout_max_time`.
        Possible values are "MAX_CALL", "MAX_OVERALL", "LAST", "AVERAGE", "ZERO".
        MAX_CALL: Use the maximum time per iteration for the current solver call as estimate.
        MAX_OVERALL: Use the maximum time per iteration over all solver calls as estimate.
        LAST: Use the time required by the last iteration as estimate.
        AVERAGE: Use an exponential moving average of the previous per iteration times as estimate (weight is currently fixed at 0.5).
        ZERO: Use 0 as estimate.
        Currently implemented for SQP only.
        Default: ZERO.
        """
        return self.__timeout_heuristic

    @timeout_heuristic.setter
    def timeout_heuristic(self, val):
        if val in ["MAX_CALL", "MAX_OVERALL", "LAST", "AVERAGE", "ZERO"]:
            self.__timeout_heuristic = val
        else:
            raise ValueError('Invalid timeout_heuristic value. Expected value in ["MAX_CALL", "MAX_OVERALL", "LAST", "AVERAGE", "ZERO"].')

    @property
    def tol(self):
        """
        NLP solver tolerance. Sets or gets the max of :py:attr:`nlp_solver_tol_eq`,
        :py:attr:`nlp_solver_tol_ineq`, :py:attr:`nlp_solver_tol_comp`
        and :py:attr:`nlp_solver_tol_stat`.
        """
        return max([self.__nlp_solver_tol_eq, self.__nlp_solver_tol_ineq,\
                    self.__nlp_solver_tol_comp, self.__nlp_solver_tol_stat])

    @tol.setter
    def tol(self, tol):
        if isinstance(tol, float) and tol > 0:
            self.__nlp_solver_tol_eq = tol
            self.__nlp_solver_tol_ineq = tol
            self.__nlp_solver_tol_stat = tol
            self.__nlp_solver_tol_comp = tol
        else:
            raise ValueError('Invalid tol value. tol must be a positive float.')

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

    @qp_tol.setter
    def qp_tol(self, qp_tol):
        if isinstance(qp_tol, float) and qp_tol > 0:
            self.__qp_solver_tol_eq = qp_tol
            self.__qp_solver_tol_ineq = qp_tol
            self.__qp_solver_tol_stat = qp_tol
            self.__qp_solver_tol_comp = qp_tol
        else:
            raise ValueError('Invalid qp_tol value. qp_tol must be a positive float.')

    @property
    def nlp_solver_tol_stat(self):
        """
        NLP solver stationarity tolerance.
        Type: float > 0
        Default: 1e-6
        """
        return self.__nlp_solver_tol_stat

    @nlp_solver_tol_stat.setter
    def nlp_solver_tol_stat(self, nlp_solver_tol_stat):
        if isinstance(nlp_solver_tol_stat, float) and nlp_solver_tol_stat > 0:
            self.__nlp_solver_tol_stat = nlp_solver_tol_stat
        else:
            raise ValueError('Invalid nlp_solver_tol_stat value. nlp_solver_tol_stat must be a positive float.')

    @property
    def nlp_solver_tol_eq(self):
        """NLP solver equality tolerance"""
        return self.__nlp_solver_tol_eq

    @nlp_solver_tol_eq.setter
    def nlp_solver_tol_eq(self, nlp_solver_tol_eq):
        if isinstance(nlp_solver_tol_eq, float) and nlp_solver_tol_eq > 0:
            self.__nlp_solver_tol_eq = nlp_solver_tol_eq
        else:
            raise ValueError('Invalid nlp_solver_tol_eq value. nlp_solver_tol_eq must be a positive float.')

    @property
    def nlp_solver_tol_min_step_norm(self):
        """
        NLP solver tolerance for minimal step norm. Solver terminates if
        step norm is below given value. If value is 0.0, then the solver does not
        test for the small step.

        Type: float
        Default: None

        If None:
        in case of FUNNEL_L1PEN_LINESEARCH: 1e-12
        otherwise: 0.0
         """
        return self.__nlp_solver_tol_min_step_norm

    @nlp_solver_tol_min_step_norm.setter
    def nlp_solver_tol_min_step_norm(self, nlp_solver_tol_min_step_norm):
        if isinstance(nlp_solver_tol_min_step_norm, float) and nlp_solver_tol_min_step_norm >= 0.0:
            self.__nlp_solver_tol_min_step_norm = nlp_solver_tol_min_step_norm
        else:
            raise ValueError('Invalid nlp_solver_tol_min_step_norm value. nlp_solver_tol_min_step_norm must be a positive float.')

    @property
    def globalization_alpha_min(self):
        """Minimal step size for globalization.

        default: None.

        If None is given:
        - in case of FUNNEL_L1PEN_LINESEARCH, value is set to 1e-17.
        - in case of MERIT_BACKTRACKING, value is set to 0.05.
        """
        return self.__globalization_alpha_min

    @globalization_alpha_min.setter
    def globalization_alpha_min(self, globalization_alpha_min):
        self.__globalization_alpha_min = globalization_alpha_min

    @property
    @deprecated(version="0.4.0", reason="Use globalization_alpha_min instead.")
    def alpha_min(self):
        """
        The option alpha_min is deprecated and has new name: globalization_alpha_min
        """
        return self.globalization_alpha_min

    @alpha_min.setter
    @deprecated(version="0.4.0", reason="Use globalization_alpha_min instead.")
    def alpha_min(self, alpha_min):
        self.globalization_alpha_min = alpha_min

    @property
    def reg_epsilon(self):
        """Epsilon for regularization, used if regularize_method in ['PROJECT', 'MIRROR', 'CONVEXIFY', 'GERSHGORIN_LEVENBERG_MARQUARDT'].

        Type: float.
        Default: 1e-4.
        """
        return self.__reg_epsilon

    @reg_epsilon.setter
    def reg_epsilon(self, reg_epsilon):
        if not isinstance(reg_epsilon, float) or reg_epsilon < 0:
            raise ValueError(f'Invalid reg_epsilon value, expected float >= 0, got {reg_epsilon}')
        self.__reg_epsilon = reg_epsilon

    @property
    def reg_max_cond_block(self):
        """Maximum condition number of each Hessian block after regularization with regularize_method in ['PROJECT', 'MIRROR'] and reg_adaptive_eps = True

        Type: float
        Default: 1e7
        """
        return self.__reg_max_cond_block

    @reg_max_cond_block.setter
    def reg_max_cond_block(self, reg_max_cond_block):
        if not isinstance(reg_max_cond_block, float) or reg_max_cond_block < 1.0:
            raise ValueError('Invalid reg_max_cond_block value, expected float > 1.0.')
        self.__reg_max_cond_block = reg_max_cond_block

    @property
    def reg_adaptive_eps(self):
        """Determines if epsilon is chosen adaptively in regularization
        used if regularize_method in ['PROJECT', 'MIRROR']

        If true, epsilon is chosen block-wise based on reg_max_cond_block.
        Otherwise, epsilon is chosen globally based on reg_epsilon.

        Type: bool
        Default: False
        """
        return self.__reg_adaptive_eps

    @reg_adaptive_eps.setter
    def reg_adaptive_eps(self, reg_adaptive_eps):
        if not isinstance(reg_adaptive_eps, bool):
            raise TypeError(f'Invalid reg_adaptive_eps value, expected bool, got {reg_adaptive_eps}')
        self.__reg_adaptive_eps = reg_adaptive_eps

    @property
    def reg_min_epsilon(self):
        """Minimum value for epsilon if regularize_method in ['PROJECT', 'MIRROR'] is used with reg_adaptive_eps.

        Type: float
        Default: 1e-8
        """
        return self.__reg_min_epsilon

    @reg_min_epsilon.setter
    def reg_min_epsilon(self, reg_min_epsilon):
        if not isinstance(reg_min_epsilon, float) or reg_min_epsilon < 0:
            raise ValueError(f'Invalid reg_min_epsilon value, expected float > 0, got {reg_min_epsilon}')
        self.__reg_min_epsilon = reg_min_epsilon

    @property
    def globalization_alpha_reduction(self):
        """Step size reduction factor for globalization MERIT_BACKTRACKING and
        FUNNEL_L1PEN_LINESEARCH

        Type: float
        Default: 0.7.
        """
        return self.__globalization_alpha_reduction

    @globalization_alpha_reduction.setter
    def globalization_alpha_reduction(self, globalization_alpha_reduction):
        self.__globalization_alpha_reduction = globalization_alpha_reduction

    @property
    @deprecated(version="0.4.0", reason="Use globalization_alpha_reduction instead.")
    def alpha_reduction(self):
        """
        The option alpha_reduction is deprecated and has new name: globalization_alpha_reduction
        """
        return self.globalization_alpha_reduction

    @alpha_reduction.setter
    @deprecated(version="0.4.0", reason="Use globalization_alpha_reduction instead.")
    def alpha_reduction(self, globalization_alpha_reduction):
        self.globalization_alpha_reduction = globalization_alpha_reduction

    @property
    def globalization_line_search_use_sufficient_descent(self):
        """
        Determines if sufficient descent (Armijo) condition is used in line search.
        Type: int; 0 or 1;
        default: 0.
        """
        return self.__globalization_line_search_use_sufficient_descent

    @globalization_line_search_use_sufficient_descent.setter
    def globalization_line_search_use_sufficient_descent(self, globalization_line_search_use_sufficient_descent):
        if globalization_line_search_use_sufficient_descent in [0, 1]:
            self.__globalization_line_search_use_sufficient_descent = globalization_line_search_use_sufficient_descent
        else:
            raise ValueError(f'Invalid value for globalization_line_search_use_sufficient_descent. Possible values are 0, 1, got {globalization_line_search_use_sufficient_descent}')

    @property
    @deprecated(version="0.4.0", reason="Use globalization_line_search_use_sufficient_descent instead.")
    def line_search_use_sufficient_descent(self):
        """
        The option line_search_use_sufficient_descent is deprecated and has new name: globalization_line_search_use_sufficient_descent
        """
        return self.globalization_line_search_use_sufficient_descent

    @line_search_use_sufficient_descent.setter
    @deprecated(version="0.4.0", reason="Use globalization_line_search_use_sufficient_descent instead.")
    def line_search_use_sufficient_descent(self, globalization_line_search_use_sufficient_descent):
        if globalization_line_search_use_sufficient_descent in [0, 1]:
            self.__globalization_line_search_use_sufficient_descent = globalization_line_search_use_sufficient_descent
        else:
            raise ValueError(f'Invalid value for globalization_line_search_use_sufficient_descent. Possible values are 0, 1, got {globalization_line_search_use_sufficient_descent}')

    @property
    def globalization_eps_sufficient_descent(self):
        """
        Factor for sufficient descent (Armijo) conditon, see also globalization_line_search_use_sufficient_descent.

        Type: float,
        Default: None.

        If None is given:
        - in case of FUNNEL_L1PEN_LINESEARCH, value is set to 1e-6.
        - in case of MERIT_BACKTRACKING, value is set to 1e-4.
        """
        return self.__globalization_eps_sufficient_descent

    @globalization_eps_sufficient_descent.setter
    def globalization_eps_sufficient_descent(self, globalization_eps_sufficient_descent):
        if isinstance(globalization_eps_sufficient_descent, float) and globalization_eps_sufficient_descent > 0:
            self.__globalization_eps_sufficient_descent = globalization_eps_sufficient_descent
        else:
            raise ValueError('Invalid globalization_eps_sufficient_descent value. globalization_eps_sufficient_descent must be a positive float.')

    @property
    @deprecated(version="0.4.0", reason="Use globalization_line_search_use_sufficient_descent instead.")
    def eps_sufficient_descent(self):
        """
        The option eps_sufficient_descent is deprecated and has new name: globalization_line_search_use_sufficient_descent
        """
        return self.globalization_line_search_use_sufficient_descent

    @eps_sufficient_descent.setter
    @deprecated(version="0.4.0", reason="Use globalization_eps_sufficient_descent instead.")
    def eps_sufficient_descent(self, globalization_eps_sufficient_descent):
        self.globalization_eps_sufficient_descent = globalization_eps_sufficient_descent

    @property
    def globalization_use_SOC(self):
        """
        Determines if second order correction (SOC) is done when using MERIT_BACKTRACKING.
        SOC is done if preliminary line search does not return full step.
        Type: int; 0 or 1;
        default: 0.
        """
        return self.__globalization_use_SOC

    @globalization_use_SOC.setter
    def globalization_use_SOC(self, globalization_use_SOC):
        if globalization_use_SOC in [0, 1]:
            self.__globalization_use_SOC = globalization_use_SOC
        else:
            raise ValueError(f'Invalid value for globalization_use_SOC. Possible values are 0, 1, got {globalization_use_SOC}')

    @property
    def globalization_full_step_dual(self):
        """
        Determines if dual variables are updated with full steps (alpha=1.0) when primal variables are updated with smaller step.

        Type: int; 0 or 1;
        default for funnel globalization: 1
        default else: 0.
        """
        return self.__globalization_full_step_dual

    @globalization_full_step_dual.setter
    def globalization_full_step_dual(self, globalization_full_step_dual):
        if globalization_full_step_dual in [0, 1]:
            self.__globalization_full_step_dual = globalization_full_step_dual
        else:
            raise ValueError(f'Invalid value for globalization_full_step_dual. Possible values are 0, 1, got {globalization_full_step_dual}')

    @property
    @deprecated(version="0.4.0", reason="Use globalization_full_step_dual instead.")
    def full_step_dual(self):
        """
        The option full_step_dual is deprecated and has new name: globalization_full_step_dual
        """
        return self.globalization_full_step_dual

    @full_step_dual.setter
    @deprecated(version="0.4.0", reason="Use globalization_full_step_dual instead.")
    def full_step_dual(self, globalization_full_step_dual):
        self.globalization_full_step_dual = globalization_full_step_dual


    @property
    def globalization_funnel_init_increase_factor(self):
        """
        Increase factor for initialization of funnel width.
        Initial funnel is max(globalization_funnel_init_upper_bound, globalization_funnel_init_increase_factor * initial_infeasibility)

        Type: float
        Default: 15.0
        """
        return self.__globalization_funnel_init_increase_factor

    @globalization_funnel_init_increase_factor.setter
    def globalization_funnel_init_increase_factor(self, globalization_funnel_init_increase_factor):
        if globalization_funnel_init_increase_factor > 1.0:
            self.__globalization_funnel_init_increase_factor = globalization_funnel_init_increase_factor
        else:
            raise ValueError(f'Invalid value for globalization_funnel_init_increase_factor. Should be > 1, got {globalization_funnel_init_increase_factor}')

    @property
    def globalization_funnel_init_upper_bound(self):
        """
        Initial upper bound for funnel width.
        Initial funnel is max(globalization_funnel_init_upper_bound, globalization_funnel_init_increase_factor * initial_infeasibility)

        Type: float
        Default: 1.0
        """
        return self.__globalization_funnel_init_upper_bound

    @globalization_funnel_init_upper_bound.setter
    def globalization_funnel_init_upper_bound(self, globalization_funnel_init_upper_bound):
        if globalization_funnel_init_upper_bound > 0.0:
            self.__globalization_funnel_init_upper_bound = globalization_funnel_init_upper_bound
        else:
             raise ValueError(f'Invalid value for globalization_funnel_init_upper_bound. Should be > 0, got {globalization_funnel_init_upper_bound}')

    @property
    def globalization_funnel_sufficient_decrease_factor(self):
        """
        Sufficient decrease factor for infeasibility in h iteration:
        trial_infeasibility <= kappa * globalization_funnel_width

        Type: float
        Default: 0.9
        """
        return self.__globalization_funnel_sufficient_decrease_factor

    @globalization_funnel_sufficient_decrease_factor.setter
    def globalization_funnel_sufficient_decrease_factor(self, globalization_funnel_sufficient_decrease_factor):
        if globalization_funnel_sufficient_decrease_factor > 0.0 and globalization_funnel_sufficient_decrease_factor < 1.0:
            self.__globalization_funnel_sufficient_decrease_factor = globalization_funnel_sufficient_decrease_factor
        else:
            raise ValueError(f'Invalid value for globalization_funnel_sufficient_decrease_factor. Should be in (0,1), got {globalization_funnel_sufficient_decrease_factor}')

    @property
    def globalization_funnel_kappa(self):
        """
        Interpolation factor for convex combination in funnel decrease function.

        Type: float
        Default: 0.9
        """
        return self.__globalization_funnel_kappa

    @globalization_funnel_kappa.setter
    def globalization_funnel_kappa(self, globalization_funnel_kappa):
        if globalization_funnel_kappa > 0.0 and globalization_funnel_kappa < 1.0:
            self.__globalization_funnel_kappa = globalization_funnel_kappa
        else:
            raise ValueError(f'Invalid value for globalization_funnel_kappa. Should be in (0,1), got {globalization_funnel_kappa}')

    @property
    def globalization_funnel_fraction_switching_condition(self):
        """
        Multiplication factor in switching condition.

        Type: float
        Default: 1e-3
        """
        return self.__globalization_funnel_fraction_switching_condition

    @globalization_funnel_fraction_switching_condition.setter
    def globalization_funnel_fraction_switching_condition(self, globalization_funnel_fraction_switching_condition):
        if globalization_funnel_fraction_switching_condition > 0.0 and globalization_funnel_fraction_switching_condition < 1.0:
            self.__globalization_funnel_fraction_switching_condition = globalization_funnel_fraction_switching_condition
        else:
            raise ValueError(f'Invalid value for globalization_funnel_fraction_switching_condition. Should be in (0,1), got {globalization_funnel_fraction_switching_condition}')

    @property
    def eval_residual_at_max_iter(self):
        """
        Determines, if the problem data is again evaluated after the last iteration
        has been performed.
        If True, the residuals are evaluated at the last iterate and convergence will be checked.
        If False, after the last iteration, the solver will terminate with max iterations
        status, although the final iterate might fulfill the convergence criteria.

        Type: bool
        Default: None

        If None is given:
        if globalization == FUNNEL_L1PEN_LINESEARCH: true
        else: false
        """
        return self.__eval_residual_at_max_iter

    @eval_residual_at_max_iter.setter
    def eval_residual_at_max_iter(self, eval_residual_at_max_iter):
        if isinstance(eval_residual_at_max_iter, bool):
            self.__eval_residual_at_max_iter = eval_residual_at_max_iter
        else:
            raise TypeError(f'Invalid datatype for eval_residual_at_max_iter. Should be bool, got {type(eval_residual_at_max_iter)}')

    @property
    def use_constraint_hessian_in_feas_qp(self):
        """
        Determines if exact/approximate Hessian of the constraints or the identity
        matrix is used as Hessian in the feasibility QP of `SQP_WITH_FEASIBLE_QP`

        Default: False
        """
        return self.__use_constraint_hessian_in_feas_qp

    @use_constraint_hessian_in_feas_qp.setter
    def use_constraint_hessian_in_feas_qp(self, use_constraint_hessian_in_feas_qp):
        if isinstance(use_constraint_hessian_in_feas_qp, bool):
            self.__use_constraint_hessian_in_feas_qp = use_constraint_hessian_in_feas_qp
        else:
            raise TypeError(f'Invalid datatype for use_constraint_hessian_in_feas_qp. Should be bool, got {type(use_constraint_hessian_in_feas_qp)}')

    @property
    def byrd_omojokon_slack_relaxation_factor(self):
        """
        Multiplication factor in Byrd-Omojokun bounds setup. Reduces ill-conditioning,
        but can allow convergence to unwanted infeasible stationary points.
        Type: double, >=1
        Default: 1.00001
        """
        return self.__byrd_omojokon_slack_relaxation_factor

    @byrd_omojokon_slack_relaxation_factor.setter
    def byrd_omojokon_slack_relaxation_factor(self, byrd_omojokon_slack_relaxation_factor):
        if isinstance(byrd_omojokon_slack_relaxation_factor, float):
            if byrd_omojokon_slack_relaxation_factor >= 1.0:
                self.__byrd_omojokon_slack_relaxation_factor = byrd_omojokon_slack_relaxation_factor
            else:
                raise ValueError(f'Invalid float for byrd_omojokon_slack_relaxation_factor. Must be >=1.0, got {byrd_omojokon_slack_relaxation_factor}')
        else:
            raise TypeError(f'Invalid datatype for search_direction_mode. Should be float, got {type(byrd_omojokon_slack_relaxation_factor)}')

    @property
    def search_direction_mode(self):
        """
        Determines how the search direction should be calculated in the initial
        iteration
        Type: string

        Possible entries are
        NOMINAL_QP, BYRD_OMOJOKUN, FEASIBILITY_QP

        Type: string
        Default: NOMINAL_QP
        Other options: BYRD_OMOJOKUN, FEASIBILITY_QP
        """
        return self.__search_direction_mode

    @search_direction_mode.setter
    def search_direction_mode(self, search_direction_mode):
        search_direction_modes = ('NOMINAL_QP', 'BYRD_OMOJOKUN', 'FEASIBILITY_QP')
        if isinstance(search_direction_mode, str):
            if search_direction_mode in search_direction_modes:
                self.__search_direction_mode = search_direction_mode
            else:
                raise ValueError(f'Invalid string for search_direction_mode. Possible search_direction_modes are'+', '.join(search_direction_modes) +  f', got {search_direction_mode}')
        else:
            raise TypeError(f'Invalid datatype for search_direction_mode. Should be str, got {type(search_direction_mode)}')

    @property
    def allow_direction_mode_switch_to_nominal(self):
        """
        Indicates if we allow switching back from BYRD_OMOJOKUN to NOMINAL_QP
        search direction mode

        Type: bool
        Default: True
        """
        return self.__allow_direction_mode_switch_to_nominal

    @allow_direction_mode_switch_to_nominal.setter
    def allow_direction_mode_switch_to_nominal(self, allow_direction_mode_switch_to_nominal):
        if isinstance(allow_direction_mode_switch_to_nominal, bool):
            self.__allow_direction_mode_switch_to_nominal = allow_direction_mode_switch_to_nominal
        else:
            raise TypeError(f'Invalid datatype for allow_direction_mode_switch_to_nominal. Should be str, got {type(allow_direction_mode_switch_to_nominal)}')

    @property
    def globalization_funnel_initial_penalty_parameter(self):
        """
        Initialization.

        Type: float
        Default: 1.0
        """
        return self.__globalization_funnel_initial_penalty_parameter

    @globalization_funnel_initial_penalty_parameter.setter
    def globalization_funnel_initial_penalty_parameter(self, globalization_funnel_initial_penalty_parameter):
        if globalization_funnel_initial_penalty_parameter >= 0.0 and globalization_funnel_initial_penalty_parameter <= 1.0:
            self.__globalization_funnel_initial_penalty_parameter = globalization_funnel_initial_penalty_parameter
        else:
            raise ValueError(f'Invalid value for globalization_funnel_initial_penalty_parameter. Should be in [0,1], got {globalization_funnel_initial_penalty_parameter}')

    @property
    def globalization_funnel_use_merit_fun_only(self):
        """
        If this options is set, the funnel globalization only checks a merit function.

        Type: bool
        Default: False
        """
        return self.__globalization_funnel_use_merit_fun_only

    @globalization_funnel_use_merit_fun_only.setter
    def globalization_funnel_use_merit_fun_only(self, globalization_funnel_use_merit_fun_only):
        if isinstance(globalization_funnel_use_merit_fun_only, bool):
            self.__globalization_funnel_use_merit_fun_only = globalization_funnel_use_merit_fun_only
        else:
            raise TypeError(f'Invalid type for globalization_funnel_use_merit_fun_only. Should be bool, got {globalization_funnel_use_merit_fun_only}')

    @property
    def nlp_solver_tol_ineq(self):
        """NLP solver inequality tolerance"""
        return self.__nlp_solver_tol_ineq

    @nlp_solver_tol_ineq.setter
    def nlp_solver_tol_ineq(self, nlp_solver_tol_ineq):
        if isinstance(nlp_solver_tol_ineq, float) and nlp_solver_tol_ineq > 0:
            self.__nlp_solver_tol_ineq = nlp_solver_tol_ineq
        else:
            raise ValueError('Invalid nlp_solver_tol_ineq value. nlp_solver_tol_ineq must be a positive float.')

    @property
    def nlp_solver_ext_qp_res(self):
        """
        Determines if residuals of QP are computed externally within NLP solver (for debugging).
        Residuals are computed on input/output of xcond-QP solver, i.e. before condensing and after expanding the QP solution.
        QP residuals are part of the statistics.
        Not supported for "SQP_WITH_FEASIBLE_QP".
        If tau_min > 0, the complementarity residuals is computed with respect to the corresponding barrier problem.

        Type: int; 0 or 1;
        Default: 0.
        """
        return self.__nlp_solver_ext_qp_res

    @nlp_solver_ext_qp_res.setter
    def nlp_solver_ext_qp_res(self, nlp_solver_ext_qp_res):
        if nlp_solver_ext_qp_res in [0, 1]:
            self.__nlp_solver_ext_qp_res = nlp_solver_ext_qp_res
        else:
            raise ValueError('Invalid nlp_solver_ext_qp_res value. nlp_solver_ext_qp_res must be in [0, 1].')

    @property
    def rti_log_residuals(self):
        """Determines if residuals are computed and logged within RTI / AS-RTI iterations (for debugging).

        Type: int; 0 or 1;
        Default: 0.
        """
        return self.__rti_log_residuals

    @rti_log_residuals.setter
    def rti_log_residuals(self, rti_log_residuals):
        if rti_log_residuals in [0, 1]:
            self.__rti_log_residuals = rti_log_residuals
        else:
            raise ValueError('Invalid rti_log_residuals value. rti_log_residuals must be in [0, 1].')

    @property
    def rti_log_only_available_residuals(self):
        """
        Relevant if rti_log_residuals is set to 1.
        If rti_log_only_available_residuals is set to 1, only residuals that do not require additional function evaluations are logged.

        Type: int; 0 or 1;
        Default: 0.
        """
        return self.__rti_log_only_available_residuals

    @rti_log_only_available_residuals.setter
    def rti_log_only_available_residuals(self, rti_log_only_available_residuals):
        if rti_log_only_available_residuals in [0, 1]:
            self.__rti_log_only_available_residuals = rti_log_only_available_residuals
        else:
            raise ValueError('Invalid rti_log_only_available_residuals value. rti_log_only_available_residuals must be in [0, 1].')

    @property
    def nlp_solver_tol_comp(self):
        """NLP solver complementarity tolerance"""
        return self.__nlp_solver_tol_comp

    @nlp_solver_tol_comp.setter
    def nlp_solver_tol_comp(self, nlp_solver_tol_comp):
        if isinstance(nlp_solver_tol_comp, float) and nlp_solver_tol_comp > 0:
            self.__nlp_solver_tol_comp = nlp_solver_tol_comp
        else:
            raise ValueError('Invalid nlp_solver_tol_comp value. nlp_solver_tol_comp must be a positive float.')

    @property
    def nlp_solver_max_iter(self):
        """
        NLP solver maximum number of iterations.
        Type: int >= 0
        Default: 100
        """
        return self.__nlp_solver_max_iter

    @nlp_solver_max_iter.setter
    def nlp_solver_max_iter(self, nlp_solver_max_iter):

        if isinstance(nlp_solver_max_iter, int) and nlp_solver_max_iter >= 0:
            self.__nlp_solver_max_iter = nlp_solver_max_iter
        else:
            raise ValueError('Invalid nlp_solver_max_iter value. nlp_solver_max_iter must be a nonnegative int.')

    @property
    def time_steps(self):
        """
        Vector of length `N_horizon` containing the time steps between the shooting nodes.
        If `None` set automatically to uniform discretization using :py:attr:`N_horizon` and :py:attr:`tf`.
        For nonuniform discretization: Either provide shooting_nodes or time_steps.
        Default: :code:`None`
        """
        return self.__time_steps

    @time_steps.setter
    def time_steps(self, time_steps):
        time_steps = check_if_nparray_and_flatten(time_steps, "time_steps")
        self.__time_steps = time_steps

    @property
    def shooting_nodes(self):
        """
        Vector of length `N_horizon + 1` containing the shooting nodes.
        If `None` set automatically to uniform discretization using :py:attr:`N_horizon` and :py:attr:`tf`.
        For nonuniform discretization: Either provide shooting_nodes or time_steps.
        Default: :code:`None`
        """
        return self.__shooting_nodes

    @shooting_nodes.setter
    def shooting_nodes(self, shooting_nodes):
        shooting_nodes = check_if_nparray_and_flatten(shooting_nodes, "shooting_nodes")
        self.__shooting_nodes = shooting_nodes

    @property
    def cost_scaling(self):
        """
        Vector with cost scaling factors of length `N_horizon` + 1.
        If `None` set automatically to [`time_steps`, 1.0].
        Default: :code:`None`
        """
        return self.__cost_scaling

    @cost_scaling.setter
    def cost_scaling(self, cost_scaling):
        cost_scaling = check_if_nparray_and_flatten(cost_scaling, "cost_scaling")
        self.__cost_scaling = cost_scaling

    @property
    def tf(self):
        """
        Prediction horizon
        Type: float > 0
        Default: :code:`None`
        """
        return self.__tf

    @tf.setter
    def tf(self, tf):
        self.__tf = tf

    @property
    def N_horizon(self):
        """
        Number of shooting intervals.
        Type: int >= 0
        Default: :code:`None`
        """
        return self.__N_horizon

    @N_horizon.setter
    def N_horizon(self, N_horizon):
        if isinstance(N_horizon, int) and N_horizon >= 0:
            self.__N_horizon = N_horizon
        else:
            raise ValueError('Invalid N_horizon value, expected non-negative integer.')

    @property
    def Tsim(self):
        """
        Time horizon for one integrator step. Automatically calculated as first time step if not provided.
        Default: :code:`None`
        """
        return self.__Tsim

    @Tsim.setter
    def Tsim(self, Tsim):
        self.__Tsim = Tsim

    @property
    def print_level(self):
        """
        Verbosity of printing.

        Type: int >= 0
        Default: 0

        Level 1: print iteration log
        Level 2: print high level debug output in funnel globalization
        Level 3: print more detailed debug output in funnel, including objective values and infeasibilities
        Level 4: print QP inputs and outputs. Please specify with max_iter how many QPs should be printed
        """
        return self.__print_level

    @print_level.setter
    def print_level(self, print_level):
        if isinstance(print_level, int) and print_level >= 0:
            self.__print_level = print_level
        else:
            raise ValueError('Invalid print_level value. print_level takes one of the values >=0.')

    @property
    def model_external_shared_lib_dir(self):
        """Path to the .so lib"""
        return self.__model_external_shared_lib_dir

    @model_external_shared_lib_dir.setter
    def model_external_shared_lib_dir(self, model_external_shared_lib_dir):
        if isinstance(model_external_shared_lib_dir, str) :
            self.__model_external_shared_lib_dir = model_external_shared_lib_dir
        else:
            raise TypeError('Invalid model_external_shared_lib_dir value. Str expected.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_dir) + '.\n\n')

    @property
    def model_external_shared_lib_name(self):
        """Name of the .so lib"""
        return self.__model_external_shared_lib_name

    @model_external_shared_lib_name.setter
    def model_external_shared_lib_name(self, model_external_shared_lib_name):
        if isinstance(model_external_shared_lib_name, str) :
            if model_external_shared_lib_name[-3:] == '.so' :
                raise ValueError('Invalid model_external_shared_lib_name value. Remove the .so extension.' \
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\n')
            else :
                self.__model_external_shared_lib_name = model_external_shared_lib_name
        else:
            raise TypeError('Invalid model_external_shared_lib_name value. Str expected.'
            + '.\n\nYou have: ' + type(model_external_shared_lib_name) + '.\n\n')

    @property
    def exact_hess_constr(self):
        """
        Used in case of hessian_approx == 'EXACT'.\n
        Can be used to turn off exact hessian contributions from the constraints module.
        """
        return self.__exact_hess_constr

    @exact_hess_constr.setter
    def exact_hess_constr(self, exact_hess_constr):
        if exact_hess_constr in [0, 1]:
            self.__exact_hess_constr = exact_hess_constr
        else:
            raise ValueError('Invalid exact_hess_constr value. exact_hess_constr takes one of the values 0, 1.')

    @property
    def exact_hess_cost(self):
        """
        Used in case of hessian_approx == 'EXACT'.\n
        Can be used to turn off exact hessian contributions from the cost module.
        """
        return self.__exact_hess_cost

    @exact_hess_cost.setter
    def exact_hess_cost(self, exact_hess_cost):
        if exact_hess_cost in [0, 1]:
            self.__exact_hess_cost = exact_hess_cost
        else:
            raise ValueError('Invalid exact_hess_cost value. exact_hess_cost takes one of the values 0, 1.')

    @property
    def exact_hess_dyn(self):
        """
        Used in case of hessian_approx == 'EXACT'.\n
        Can be used to turn off exact hessian contributions from the dynamics module.
        """
        return self.__exact_hess_dyn

    @exact_hess_dyn.setter
    def exact_hess_dyn(self, exact_hess_dyn):
        if exact_hess_dyn in [0, 1]:
            self.__exact_hess_dyn = exact_hess_dyn
        else:
            raise ValueError('Invalid exact_hess_dyn value. exact_hess_dyn takes one of the values 0, 1.')

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

    @fixed_hess.setter
    def fixed_hess(self, fixed_hess):
        if fixed_hess in [0, 1]:
            self.__fixed_hess = fixed_hess
        else:
            raise ValueError('Invalid fixed_hess value. fixed_hess takes one of the values 0, 1.')

    @property
    def ext_cost_num_hess(self):
        """
        Determines if custom hessian approximation for cost contribution is used (> 0).\n
        Or if hessian contribution is evaluated exactly using CasADi external function (=0 - default).
        """
        return self.__ext_cost_num_hess

    @ext_cost_num_hess.setter
    def ext_cost_num_hess(self, ext_cost_num_hess):
        if ext_cost_num_hess in [0, 1]:
            self.__ext_cost_num_hess = ext_cost_num_hess
        else:
            raise ValueError('Invalid ext_cost_num_hess value. ext_cost_num_hess takes one of the values 0, 1.')

    @property
    def cost_discretization(self):
        """
        Cost discretization: string in {'EULER', 'INTEGRATOR'}.
        Default: 'EULER'
        'EULER': cost is evaluated at shooting nodes
        'INTEGRATOR': cost is integrated over the shooting intervals - only supported for IRK integrator
        """
        return self.__cost_discretization

    @cost_discretization.setter
    def cost_discretization(self, cost_discretization):
        if cost_discretization in COST_DISCRETIZATION_TYPES:
            self.__cost_discretization = cost_discretization
        else:
            raise ValueError('Invalid cost_discretization value. Possible values are:\n\n' \
                    + ',\n'.join(COST_DISCRETIZATION_TYPES) + '.\n\nYou have: ' + cost_discretization + '.')

    @property
    def with_solution_sens_wrt_params(self):
        """
        Flag indicating whether solution sensitivities wrt. parameters can be computed.
        """
        return self.__with_solution_sens_wrt_params

    @with_solution_sens_wrt_params.setter
    def with_solution_sens_wrt_params(self, with_solution_sens_wrt_params):
        if isinstance(with_solution_sens_wrt_params, bool):
            self.__with_solution_sens_wrt_params = with_solution_sens_wrt_params
        else:
            raise TypeError('Invalid with_solution_sens_wrt_params value. Expected bool.')

    @property
    def with_value_sens_wrt_params(self):
        """
        Flag indicating whether value function sensitivities wrt. parameters can be computed.
        """
        return self.__with_value_sens_wrt_params

    @with_value_sens_wrt_params.setter
    def with_value_sens_wrt_params(self, with_value_sens_wrt_params):
        if isinstance(with_value_sens_wrt_params, bool):
            self.__with_value_sens_wrt_params = with_value_sens_wrt_params
        else:
            raise TypeError('Invalid with_value_sens_wrt_params value. Expected bool.')

    @property
    @deprecated(version="0.4.0", reason="Use with_batch_functionality instead and pass the number of threads directly to the BatchSolver.")
    def num_threads_in_batch_solve(self):
        """
        Integer indicating how many threads should be used within the batch solve.
        If more than one thread should be used, the solver is compiled with openmp.
        Default: 1.
        """
        return self.__num_threads_in_batch_solve

    @num_threads_in_batch_solve.setter
    @deprecated(version="0.4.0", reason="Set the flag with_batch_functionality instead and pass the number of threads directly to the BatchSolver.")
    def num_threads_in_batch_solve(self, num_threads_in_batch_solve):
        if isinstance(num_threads_in_batch_solve, int) and num_threads_in_batch_solve > 0:
            self.__num_threads_in_batch_solve = num_threads_in_batch_solve
        else:
            raise ValueError('Invalid num_threads_in_batch_solve value. num_threads_in_batch_solve must be a positive integer.')

    @property
    def with_batch_functionality(self):
        """
        Whether the AcadosOcpBatchSolver can be used.
        In this case, the solver is compiled with openmp.
        Default: False.
        """
        return self.__with_batch_functionality


    @with_batch_functionality.setter
    def with_batch_functionality(self, with_batch_functionality):
        if isinstance(with_batch_functionality, bool):
            self.__with_batch_functionality = with_batch_functionality
        else:
            raise TypeError('Invalid with_batch_functionality value. Expected bool.')


    def set(self, attr, value):
        setattr(self, attr, value)

    def _ensure_solution_sensitivities_available(self, parametric: bool = True, has_custom_hess: bool = False):
        if not self.qp_solver in ['FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_HPIPM']:
            raise NotImplementedError("Parametric sensitivities are only available with HPIPM as QP solver.")

        if not (
            self.hessian_approx == 'EXACT' and
            self.regularize_method == 'NO_REGULARIZE' and
            self.levenberg_marquardt == 0 and
            self.exact_hess_constr == 1 and
            self.exact_hess_cost == 1 and
            self.exact_hess_dyn == 1 and
            self.fixed_hess == 0 and
            has_custom_hess is False
        ):
            raise ValueError("Parametric sensitivities are only correct if an exact Hessian is used!")

        if parametric and not self.with_solution_sens_wrt_params:
            raise ValueError("Parametric sensitivities are only available if with_solution_sens_wrt_params is set to True.")

        if self.qpscaling_scale_constraints != "NO_CONSTRAINT_SCALING" or self.qpscaling_scale_objective != "NO_OBJECTIVE_SCALING":
            raise ValueError("Parametric sensitivities are only available if no scaling is applied to the QP.")

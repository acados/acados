class ocp_nlp_dims:
    def __init__(self):
        self._nx = None  # number of states
        self._nu = None  # number of inputs
        self._N  = None  # prediction horizon 

    @property
    def nx(self):
        return self._nx

    @property
    def nu(self):
        return self._nu

    @property
    def N(self):
        return self._N

    @nx.setter
    def nx(self, nx):
        if type(nx) == int and nx > 0:
            self._nx = nx
        else:
            raise Exception('Invalid nx value. Exiting.')

    @nu.setter
    def nu(self, nu):
        if type(nu) == int and nu > 0:
            self._nu = nu
        else:
            raise Exception('Invalid nu value. Exiting.')

    @N.setter
    def N(self, N):
        if type(N) == int and N > 0:
            self._N = N
        else:
            raise Exception('Invalid N value. Exiting.')

class ocp_nlp_solver_config:
    def __init__(self):
        self._qp_solver      = 'PARTIAL_CONDENSING_HPIPM' # qp solver to be used in the NLP solver
        self._hessian_approx = 'GAUSS_NEWTON' # hessian approximation
        self._integrator_type = 'ERK' # integrator type

    @property
    def qp_solver(self):
        return self._qp_solver

    @property
    def hessian_approx(self):
        return self._hessian_approx

    @property
    def integrator_type(self):
        return self._integrator_type

    @qp_solver.setter
    def qp_solver(self, qp_solver):
        qp_solvers = ('PARTIAL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_QPOASES', \
                'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_QPOASES')

        if type(qp_solver) == str and qp_solver in qp_solvers:
            self._qp_solver = qp_solver
        else:
            raise Exception('Invalid qp_solver value. Possible values are:\n\n' \
                    + ',\n'.join(qp_solvers) + '.\n\nYou have: ' + qp_solver + '.\n\nExiting.')

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
        integrator_types = ('ERK')

        if type(integrator_type) == str and integrator_type in integrator_types:
            self._integrator_type = integrator_type
        else:
            raise Exception('Invalid integrator_type value. Possible values are:\n\n' \
                    + ',\n'.join(integrator_types) + '.\n\nYou have: ' + integrator_type + '.\n\nExiting.')
class ocp_nlp_constant:
    def __init__(self):
        self.name  = None # constant name
        self.value = None # constant value

class ocp_nlp_render_arguments:
    def __init__(self):
        self.dims = ocp_nlp_dims()
        self.solver_config = ocp_nlp_solver_config()
        self.model_name = None 
        self.constants = []
        self.acados_include_path = []

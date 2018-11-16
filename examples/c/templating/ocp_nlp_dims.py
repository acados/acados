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
        self.qp_solver      = None # qp solve to be used in the NLP solver
        self.hessian_approx = None # hessian approximation

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

class ocp_nlp_dims:
    def __init__(self):

        self.nx = None  # number of states
        self.nu = None  # number of inputs
        self.N  = None  # prediction horizon 

class ocp_nlp_solver_config:
    def __init__(self):
        self.qp_solver = None # qp solve to be used in the NLP solver

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

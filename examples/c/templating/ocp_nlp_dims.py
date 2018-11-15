class ocp_nlp_dims:
    def __init__(self):

        self.nx = None  # number of states
        self.nu = None  # number of inputs
        self.ni = None  # number of two-sided inequality constraints TODO make one-sided ???
        self.N  = None  # prediction horizon 

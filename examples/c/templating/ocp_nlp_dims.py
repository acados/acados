class ocp_nlp_dims:
    def __init__(self):
        # void **cost;
        # void **dynamics;
        # void **constraints;
        # ocp_qp_dims *qp_solver;  // xcond solver instead ??

        self.nv = None  # number of primal variables (states+controls+slacks)
        self.nx = None  # number of states
        self.nu = None  # number of inputs
        self.ni = None  # number of two-sided inequality constraints TODO make one-sided ???
        self.nz = None  # number of algebraic variables 
        self.N  = None  # prediction horizon 

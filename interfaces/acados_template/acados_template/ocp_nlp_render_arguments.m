classdef ocp_nlp_dims < handle
    properties
        nx   % number of states
        nz   % number of algebraic variables
        nu   % number of inputs
        np   % number of parameters
        ny   % number of residuals in Lagrange term
        nyN  % number of residuals in Mayer term
        npd  % number of positive definite constraints
        npdN % number of positive definite constraints in last stage
        nh   % number of nonlinear constraints
        nhN  % number of nonlinear constraints in last stage
        nbx  % number of state bounds 
        nbu  % number of input bounds
        ng   % number of general constraints
        nbxN % number of state bounds in last stage 
        ngN  % number of general constraints in last stage
        N    % prediction horizon 
    end
    methods 
        function obj = ocp_nlp_dims()
            obj.nx   = []; 
            obj.nz   = 0; 
            obj.nu   = []; 
            obj.np   = 0
            obj.ny   = []; 
            obj.nyN  = []; 
            obj.npd  = 0;
            obj.npdN = 0; 
            obj.nh   = 0;
            obj.nhN  = 0;
            obj.nbx  = 0;
            obj.nbu  = 0;
            obj.ng   = 0;
            obj.nbxN = 0;
            obj.ngN  = 0;
            obj.N    = [];
        end
    end
end

classdef ocp_nlp_cost < handle
    % linear least-squares cost: || Vx*x + Vu*x + Vz*z ||^2_W
    properties
        % Lagrange term
        W      % weight matrix
        Vx     % x matrix coefficient
        Vu     % u matrix coefficient
        Vz     % z matrix coefficient
        yref   % reference
        % Mayer term
        WN     % weight matrix
        VxN    % x matrix coefficient
        yrefN  % reference
    end
    methods
        function obj = ocp_nlp_cost()
            obj.W     = [];  
            obj.Vx    = [];
            obj.Vu    = [];
            obj.Vz    = [];
            obj.yref  = [];
            obj.WN    = [];
            obj.VxN   = [];
            obj.yrefN = [];
        end
    end
end

classdef ocp_nlp_constraints < handle
    properties
        lbx     % lower bounds on x
        lbu     % lower bounds on u
        idxbx   % indexes of bounds on x 
        ubx     % upper bounds on x 
        ubu     % upper bounds on u 
        idxbu   % indexes of bounds on u
        lg      % lower bound for general inequalities 
        ug      % upper bound for general inequalities 
        lh      % lower bound for nonlinear inequalities 
        uh      % upper bound for nonlinear inequalities 
        D       % D matrix in lg <= D * u + C * x <= ug
        C       % C matrix in lg <= D * u + C * x <= ug
        lbxN    % lower bounds on x at t=T 
        ubxN    % upper bounds on x at t=T 
        idxbxN  % indexes for bounds on x at t=T 
        CN      % C matrix at t=T 
        lgN     % lower bound on general inequalities at t=T 
        ugN     % upper bound on general inequalities at t=T 
        lhN     % lower bound on nonlinear inequalities at t=T 
        uhN     % upper bound on nonlinear inequalities at t=T 
        x0      % initial state 
    end
    methods
        function obj = ocp_nlp_constraints()
            obj.lbx    = [];  
            obj.lbu    = [];
            obj.idxbx  = [];
            obj.ubx    = [];
            obj.ubu    = [];
            obj.idxbu  = [];
            obj.lg     = [];
            obj.ug     = [];
            obj.lh     = [];
            obj.uh     = [];
            obj.D      = [];
            obj.C      = [];
            obj.lbxN   = [];
            obj.ubxN   = [];
            obj.idxbxN = [];
            obj.CN     = [];
            obj.lgN    = [];
            obj.ugN    = [];
            obj.lhN    = [];
            obj.uhN    = [];
            obj.x0     = [];
        end
    end
end

classdef ocp_nlp_solver_config < handle
    properties
        qp_solver        %  qp solver to be used in the NLP solver
        hessian_approx   %  hessian approximation
        integrator_type  %  integrator type
        tf               %  prediction horizon
        nlp_solver_type  %  NLP solver 
    end
    methods
        function obj = ocp_nlp_solver_config()
            obj.qp_solver       = "PARTIAL_CONDENSING_HPIPM"; 
            obj.hessian_approx  = "GAUSS_NEWTON"; 
            obj.integrator_type = "ERK"; 
            obj.tf              = []; 
            obj.nlp_solver_type = "SQP_RTI"; 
        end
    end
end

classdef ocp_nlp_constant < handle
    properties
        name  % constant name
        value % constant value
    end
    methods
        function obj = ocp_nlp_constant()
            obj.name  = [];
            obj.value = [];    
        end
    end
end

classdef ocp_nlp_render_arguments:
    properties
        dims 
        cost 
        constraints 
        solver_config 
        model_name 
        con_p_name 
        con_pN_name 
        con_h_name 
        con_hN_name 
        constants 
        acados_include_path 
        acados_lib_path 
    end
    methods 
        function obj = ocp_nlp_render_arguments()
            obj.dims = ocp_nlp_dims(); 
            obj.cost = ocp_nlp_cost();
            obj.constraints = ocp_nlp_constraints();
            obj.solver_config = ocp_nlp_solver_config(); 
            obj.model_name = []; 
            obj.con_p_name = [];
            obj.con_pN_name = [];
            obj.con_h_name = [];
            obj.con_hN_name = [];
            obj.constants = [];
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
        end
    end
end


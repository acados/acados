classdef ocp_nlp_dims < handle
    properties
        nx    % number of states
        nz    % number of algebraic variables
        nu    % number of inputs
        np    % number of parameters
        ny    % number of residuals in Lagrange term
        nyN   % number of residuals in Mayer term
        npd   % number of positive definite constraints
        npdN  % number of positive definite constraints at t=T
        nh    % number of nonlinear constraints
        nhN   % number of nonlinear constraints at t=T
        nbx   % number of state bounds 
        nbxN  % number of state bounds at t=T 
        nbu   % number of input bounds
        nsbx  % number of soft state bounds 
        nsbu  % number of soft state bounds 
        nsbxN % number of state bounds at t=T 
        ns    % total number of soft bounds 
        nsN   % total number of soft bounds at t=T 
        ng    % number of general constraints
        ngN   % number of general constraints at t=T
        N     % prediction horizon 
    end
    methods 
        function obj = ocp_nlp_dims()
            obj.nx    = []; 
            obj.nz    = 0; 
            obj.nu    = []; 
            obj.np    = 0;
            obj.ny    = []; 
            obj.nyN   = []; 
            obj.npd   = 0;
            obj.npdN  = 0; 
            obj.nh    = 0;
            obj.nhN   = 0;
            obj.nbx   = 0;
            obj.nbu   = 0;
            obj.nbxN  = 0;
            obj.nsbx  = 0;
            obj.nsbu  = 0;
            obj.nsbxN = 0;
            obj.ns    = 0;
            obj.nsN   = 0;
            obj.ng    = 0;
            obj.ngN   = 0;
            obj.N     = [];
        end
    end
end


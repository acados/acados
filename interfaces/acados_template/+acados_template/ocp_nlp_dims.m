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
            obj.np   = 0;
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


classdef ocp_nlp_dims_json < handle
    properties
        nx     % number of states
        nz     % number of algebraic variables
        nu     % number of inputs
        np     % number of parameters
        ny     % number of residuals in Lagrange term
        ny_e   % number of residuals in Mayer term
        npd    % number of positive definite constraints
        npd_e  % number of positive definite constraints at t=T
        nh     % number of nonlinear constraints
        nh_e   % number of nonlinear constraints at t=T
        nbx    % number of state bounds 
        nbx_e  % number of state bounds at t=T 
        nbu    % number of input bounds
        nsbx   % number of soft state bounds 
        nsbu   % number of soft state bounds 
        nsbx_e % number of state bounds at t=T 
        ns     % total number of soft bounds 
        ns_e   % total number of soft bounds at t=T 
        nsh    % number of soft bounds on nonlinear constraints 
        nsh_e  % number of soft bounds on nonlinear constraints at t=T 
        ng     % number of general constraints
        ng_e   % number of general constraints at t=T
        N      % prediction horizon 
        % Declare convex over nonlinear stuff that should be implemented in MEX
        % TODO..
        nphi   % dimension of convex outer part for convex over nonlinear constraint (BGP)
        nphi_e % dimension of convex outer part for convex over nonlinear constraint (BGP) at t=T
        nsphi  % number of softend convex over nonlinear constraints (BGP)
        nsphi_e % number of softend convex over nonlinear constraints (BGP) at t=T
        nr %
        nr_e %

    end
    methods 
        function obj = ocp_nlp_dims_json()
            obj.nx    = []; 
            obj.nz    = 0; 
            obj.nu    = []; 
            obj.np    = 0;
            obj.ny    = []; 
            obj.ny_e   = []; 
            obj.npd   = 0;
            obj.npd_e  = 0; 
            obj.nh    = 0;
            obj.nh_e   = 0;
            obj.nbx   = 0;
            obj.nbu   = 0;
            obj.nbx_e  = 0;
            obj.nsbx  = 0;
            obj.nsbu  = 0;
            obj.nsbx_e = 0;
            obj.ns    = 0;
            obj.ns_e   = 0;
            obj.nsh   = 0;
            obj.nsh_e  = 0;
            obj.ng    = 0;
            obj.ng_e   = 0;
            obj.nphi = 0;
            obj.nsphi = 0;
            obj.nphi_e = 0;
            obj.nsphi_e = 0;
            obj.nr = 0;
            obj.nr_e = 0;
            obj.N     = [];
        end
    end
end


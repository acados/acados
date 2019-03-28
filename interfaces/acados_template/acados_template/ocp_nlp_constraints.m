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


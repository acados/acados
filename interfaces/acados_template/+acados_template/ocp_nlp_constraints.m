classdef ocp_nlp_constraints < handle
    properties
        % bounds on x and u
        lbx     % lower bounds on x
        lbu     % lower bounds on u
        idxbx   % indexes of bounds on x 
        ubx     % upper bounds on x 
        ubu     % upper bounds on u 
        idxbu   % indexes of bounds on u
        % bounds on x at t=T
        lbxN    % lower bounds on x at t=T 
        ubxN    % upper bounds on x at t=T 
        idxbxN  % indexes for bounds on x at t=T 
        % soft bounds on x and u
        lsbx    % soft lower bounds on x
        lsbu    % soft lower bounds on u
        usbx    % soft upper bounds on x 
        usbu    % soft upper bounds on u 
        idxsbx  % indexes of soft bounds on x 
        idxsbu  % indexes of soft bounds on u
        % soft bounds on x and u at t=T
        lsbxN   % soft lower bounds on x at t=T
        usbxN   % soft upper bounds on x at t=T
        idxsbxN % indexes of soft bounds on x at t=T 
        % polytopic constraints 
        D       % D matrix in lg <= D * u + C * x <= ug
        C       % C matrix in lg <= D * u + C * x <= ug
        lg      % lower bound for general inequalities 
        ug      % upper bound for general inequalities 
        % polytopic constraints at t=T 
        CN      % C matrix at t=T 
        lgN     % lower bound on general inequalities at t=T 
        ugN     % upper bound on general inequalities at t=T 
        % nonlinear constraints
        lh      % lower bound for nonlinear inequalities 
        uh      % upper bound for nonlinear inequalities 
        % nonlinear constraints at t=T
        lhN     % lower bound on nonlinear inequalities at t=T 
        uhN     % upper bound on nonlinear inequalities at t=T 
        x0      % initial state 
    end
    methods
        function obj = ocp_nlp_constraints()
            obj.lbx     = [];  
            obj.lbu     = [];
            obj.idxbx   = [];
            obj.ubx     = [];
            obj.ubu     = [];
            obj.idxbu   = [];
            obj.lsbx    = [];  
            obj.lsbu    = [];
            obj.idxsbx  = [];
            obj.usbx    = [];
            obj.usbu    = [];
            obj.idxsbu  = [];
            obj.lsbxN   = [];  
            obj.lsbuN   = [];
            obj.idxsbxN = [];
            obj.lg      = [];
            obj.ug      = [];
            obj.lh      = [];
            obj.uh      = [];
            obj.D       = [];
            obj.C       = [];
            obj.lbxN    = [];
            obj.ubxN    = [];
            obj.idxbxN  = [];
            obj.CN      = [];
            obj.lgN     = [];
            obj.ugN     = [];
            obj.lhN     = [];
            obj.uhN     = [];
            obj.x0      = [];
        end
    end
end


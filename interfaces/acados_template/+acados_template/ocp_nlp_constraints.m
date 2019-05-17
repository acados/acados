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
        lbx_e    % lower bounds on x at t=T 
        ubx_e    % upper bounds on x at t=T 
        idxbx_e  % indexes for bounds on x at t=T 
        % soft bounds on x and u
        lsbx    % soft lower bounds on x
        lsbu    % soft lower bounds on u
        usbx    % soft upper bounds on x 
        usbu    % soft upper bounds on u 
        idxsbx  % indexes of soft bounds on x 
        idxsbu  % indexes of soft bounds on u
        % soft bounds on x and u at t=T
        lsbx_e   % soft lower bounds on x at t=T
        usbx_e   % soft upper bounds on x at t=T
        idxsbx_e % indexes of soft bounds on x at t=T 
        % polytopic constraints 
        D       % D matrix in lg <= D * u + C * x <= ug
        C       % C matrix in lg <= D * u + C * x <= ug
        lg      % lower bound for general inequalities 
        ug      % upper bound for general inequalities 
        % polytopic constraints at t=T 
        C_e      % C matrix at t=T 
        lg_e     % lower bound on general inequalities at t=T 
        ug_e     % upper bound on general inequalities at t=T 
        % nonlinear constraints
        lh      % lower bound for nonlinear inequalities 
        uh      % upper bound for nonlinear inequalities 
        % nonlinear constraints at t=T
        lh_e     % lower bound on nonlinear inequalities at t=T 
        uh_e     % upper bound on nonlinear inequalities at t=T 
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
            obj.lsbx_e   = [];  
            obj.idxsbx_e = [];
            obj.lg      = [];
            obj.ug      = [];
            obj.lh      = [];
            obj.uh      = [];
            obj.D       = [];
            obj.C       = [];
            obj.lbx_e    = [];
            obj.ubx_e    = [];
            obj.idxbx_e  = [];
            obj.C_e      = [];
            obj.lg_e     = [];
            obj.ug_e     = [];
            obj.lh_e     = [];
            obj.uh_e     = [];
            obj.x0      = [];
        end
    end
end


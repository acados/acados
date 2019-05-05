classdef ocp_nlp_cost < handle
    % linear least-squares cost: || Vx*x + Vu*x + Vz*z ||^2_W
    properties
        % Lagrange term
        W      % weight matrix
        Vx     % x matrix coefficient
        Vu     % u matrix coefficient
        Vz     % z matrix coefficient
        yref   % reference
        Zl     % Hessian wrt lower slack 
        Zu     % Hessian wrt upper slack 
        zl     % gradient wrt lower slack 
        zu     % gradient wrt upper slack 
        % Mayer term
        WN     % weight matrix
        VxN    % x matrix coefficient
        yrefN  % reference
        ZlN    % Hessian wrt lower slack 
        ZuN    % Hessian wrt upper slack 
        zlN    % gradient wrt lower slack 
        zuN    % gradient wrt upper slack 
    end
    methods
        function obj = ocp_nlp_cost()
            obj.W     = [];  
            obj.Vx    = [];
            obj.Vu    = [];
            obj.Vz    = [];
            obj.yref  = [];
            obj.Zl    = [];
            obj.Zu    = [];
            obj.zl    = [];
            obj.zu    = [];
            obj.WN    = [];
            obj.VxN   = [];
            obj.yrefN = [];
            obj.ZlN   = [];
            obj.ZuN   = [];
            obj.zlN   = [];
            obj.zuN   = [];
        end
    end
end


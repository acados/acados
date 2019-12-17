classdef ocp_nlp_cost_json < handle
    % linear least-squares cost: || Vx*x + Vu*x + Vz*z ||^2_W
    properties
        % Lagrange term
        cost_type   % cost type
        W           % weight matrix
        Vx          % x matrix coefficient
        Vu          % u matrix coefficient
        Vz          % z matrix coefficient
        yref        % reference
        Zl          % Hessian wrt lower slack 
        Zu          % Hessian wrt upper slack 
        zl          % gradient wrt lower slack 
        zu          % gradient wrt upper slack 
        % Mayer term
        cost_type_e % cost type
        W_e         % weight matrix
        Vx_e        % x matrix coefficient
        yref_e      % reference
        Zl_e        % Hessian wrt lower slack 
        Zu_e        % Hessian wrt upper slack 
        zl_e        % gradient wrt lower slack 
        zu_e        % gradient wrt upper slack 
    end
    methods
        function obj = ocp_nlp_cost()
            obj.cost_type   = 'LINEAR_LS';  
            obj.W           = [];  
            obj.Vx          = [];
            obj.Vu          = [];
            obj.Vz          = [];
            obj.yref        = [];
            obj.Zl          = [];
            obj.Zu          = [];
            obj.zl          = [];
            obj.zu          = [];
            obj.cost_type_e = 'LINEAR_LS';  
            obj.W_e         = [];
            obj.Vx_e        = [];
            obj.yref_e      = [];
            obj.Zl_e        = [];
            obj.Zu_e        = [];
            obj.zl_e        = [];
            obj.zu_e        = [];
        end
    end
end


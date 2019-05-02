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


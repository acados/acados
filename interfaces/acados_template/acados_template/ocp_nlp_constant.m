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


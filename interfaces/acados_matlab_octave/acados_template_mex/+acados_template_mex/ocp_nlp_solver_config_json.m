classdef ocp_nlp_solver_config_json < handle
    properties
        qp_solver        %  qp solver to be used in the NLP solver
        hessian_approx   %  hessian approximation
        integrator_type  %  integrator type
        tf               %  prediction horizon
        nlp_solver_type  %  NLP solver 
    end
    methods
        function obj = ocp_nlp_solver_config_json()
            obj.qp_solver       = 'PARTIAL_CONDENSING_HPIPM';
            obj.hessian_approx  = 'GAUSS_NEWTON';
            obj.integrator_type = 'ERK';
            obj.tf              = [];
            obj.nlp_solver_type = 'SQP_RTI';
        end
    end
end

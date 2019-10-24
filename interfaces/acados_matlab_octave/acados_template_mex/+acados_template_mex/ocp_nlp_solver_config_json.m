classdef ocp_nlp_solver_config_json < handle
    properties
        qp_solver              %  qp solver to be used in the NLP solver
        hessian_approx         %  hessian approximation
        integrator_type        %  integrator type
        tf                     %  prediction horizon
        nlp_solver_type        %  NLP solver
        sim_method_num_steps   %  number of steps in integrator
        sim_method_num_stages  %  size of butcher tableau
    end
    methods
        function obj = ocp_nlp_solver_config_json()
            obj.qp_solver       = 'PARTIAL_CONDENSING_HPIPM';
            obj.hessian_approx  = 'GAUSS_NEWTON';
            obj.integrator_type = 'ERK';
            obj.tf              = [];
            obj.nlp_solver_type = 'SQP_RTI';
            obj.sim_method_num_steps = 1;
            obj.sim_method_num_stages = 2;
        end
    end
end

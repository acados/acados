classdef acados_ocp_nlp_json < handle
    properties
        dims 
        cost 
        constraints 
        solver_config 
        model_name 
        con_p_name 
        con_p_e_name 
        con_h_name 
        con_h_e_name 
        constants 
        acados_include_path 
        acados_lib_path 
    end
    methods 
        function obj = acados_ocp_nlp_json()
            obj.dims = acados_template_mex.ocp_nlp_dims_json(); 
            obj.cost = acados_template_mex.ocp_nlp_cost_json();
            obj.constraints = acados_template_mex.ocp_nlp_constraints_json();
            obj.solver_config = acados_template_mex.ocp_nlp_solver_config_json(); 
            obj.model_name = []; 
            obj.con_p_name = [];
            obj.con_p_e_name = [];
            obj.con_h_name = [];
            obj.con_h_e_name = [];
            obj.constants = [];
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
        end
    end
end


classdef acados_ocp_nlp < handle
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
        function obj = acados_ocp_nlp()
            obj.dims = acados_template.ocp_nlp_dims(); 
            obj.cost = acados_template.ocp_nlp_cost();
            obj.constraints = acados_template.ocp_nlp_constraints();
            obj.solver_config = acados_template.ocp_nlp_solver_config(); 
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


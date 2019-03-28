classdef ocp_nlp_render_arguments < handle
    properties
        dims 
        cost 
        constraints 
        solver_config 
        model_name 
        con_p_name 
        con_pN_name 
        con_h_name 
        con_hN_name 
        constants 
        acados_include_path 
        acados_lib_path 
    end
    methods 
        function obj = ocp_nlp_render_arguments()
            obj.dims = ocp_nlp_dims(); 
            obj.cost = ocp_nlp_cost();
            obj.constraints = ocp_nlp_constraints();
            obj.solver_config = ocp_nlp_solver_config(); 
            obj.model_name = []; 
            obj.con_p_name = [];
            obj.con_pN_name = [];
            obj.con_h_name = [];
            obj.con_hN_name = [];
            obj.constants = [];
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
        end
    end
end


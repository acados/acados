classdef acados_ocp_nlp_json < handle
    properties
        dims 
        cost 
        constraints 
        solver_config 
        model 
        con_p 
        con_p_e 
        con_h 
        con_h_e 
        acados_include_path 
        acados_lib_path 
    end
    methods 
        function obj = acados_ocp_nlp_json()
            obj.dims = acados_template_mex.ocp_nlp_dims_json(); 
            obj.cost = acados_template_mex.ocp_nlp_cost_json();
            obj.constraints = acados_template_mex.ocp_nlp_constraints_json();
            obj.solver_config = acados_template_mex.ocp_nlp_solver_config_json(); 
            obj.model = acados_template_mex.acados_dae(); 
            obj.con_p = acados_template_mex.acados_constraint();
            obj.con_p_e = acados_template_mex.acados_constraint();
            obj.con_h = acados_template_mex.acados_constraint();
            obj.con_h_e = acados_template_mex.acados_constraint();
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
        end
    end
end


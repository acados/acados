function simulink_opts = add_sparse_param_port_simulink(simulink_opts, idx_p, port_name, stage_idx_0, stage_idx_e)
    %% add_sparse_param_port_simulink:
    % allows one to specify information for an input port of the simulink block corresponding to an acados OCP solver.
    %
    % inputs:
    % idx_p is a 0-based vector of parameter indices to be updated by the port.
    % port_name is used to identify the port and print information
    % stage_idx_0 is the first stage for which the parameters should be updated by the port
    % stage_idx_e is the last stage for which the parameters should be updated by the port

    % sanity checks:
    if stage_idx_0 > stage_idx_e
        error("stage_idx_0 > stage_idx_e")
    end

    input_name = strcat('sparse_parameter_', port_name);
    if ~isfield(simulink_opts, 'customizable_inputs')
        simulink_opts.customizable_inputs = struct();
    end
    simulink_opts.customizable_inputs = setfield(simulink_opts.customizable_inputs, input_name, ...
            struct('parameter_indices', idx_p, 'stage_idx_0', stage_idx_0, 'stage_idx_e', stage_idx_e));

    % NOTE: putting this logic before somehow does not work in Maltab...
    if length(idx_p) == 1
        simulink_opts.customizable_inputs.(input_name).parameter_indices = reshape(num2cell(idx_p), [1, 1]);
    end
end

function simulink_opts = add_sparse_param_port_simulink(simulink_opts, idx_p, port_name, stage_idx_0, stage_idx_e)
    %% add_sparse_param_port_simulink:
    % DEPRECATED: use simulink_opts.add_sparse_param_port(idx_p, port_name, stage_idx_0, stage_idx_e) instead.
    % This function will be removed in a future release.
    warning(['add_sparse_param_port_simulink() is deprecated and will be removed in a future release. ', ...
        'Use simulink_opts.add_sparse_param_port(idx_p, port_name, stage_idx_0, stage_idx_e) instead.']);

    if ~isa(simulink_opts, 'AcadosOcpSimulinkOptions')
        error('simulink_opts must be an AcadosOcpSimulinkOptions object.');
    end

    simulink_opts.add_sparse_param_port(idx_p, port_name, stage_idx_0, stage_idx_e);
end

classdef ZoroDescription < handle
    properties

        backoff_scaling_gamma = 1.0

        feedback_optimization_mode = 'CONSTANT_FEEDBACK'

        parameter_uncertainty_mode = 'CONSTANT'
        % Type of parameter uncertainty propagation used in zoRO covariance recursion.
        %
        % String in: 'CONSTANT', 'IID', 'NONE'
        % - CONSTANT: fixed parameter error over the horizon (Π Σ_p Π^T formulation)
        % - IID: stepwise (i.i.d.) parameter uncertainty (S_p Σ_p S_p^T per stage)
        % - NONE: ignore parameter uncertainty

        fdbk_K_mat = []

        riccati_Q_const = []
        riccati_Q_const_e = []
        riccati_R_const = []
        riccati_S_const = []

        unc_jac_G_mat = []
        P0_mat = []
        W_mat = []
        Sigma_p_mat = []


        idx_lbx_t = []
        idx_ubx_t = []
        idx_lbx_e_t = []
        idx_ubx_e_t = []
        idx_lbu_t = []
        idx_ubu_t = []
        idx_lg_t = []
        idx_ug_t = []
        idx_lg_e_t = []
        idx_ug_e_t = []
        idx_lh_t = []
        idx_uh_t = []
        idx_lh_e_t = []
        idx_uh_e_t = []




        % streaming options
        input_P0_diag      = false
        input_P0           = true
        input_W_diag       = false
        input_W_add_diag   = false
        input_Sigma_p_diag = false
        input_Sigma_p      = false

        output_P_matrices = false
        output_riccati_t = false

    % properties (Access = private)
    % kind of private, but need to be dumped to json
        np
        nw
        nlbx_t
        nubx_t
        nlbx_e_t
        nubx_e_t
        nlbu_t
        nubu_t
        nlg_t
        nug_t
        nlg_e_t
        nug_e_t
        nlh_t
        nuh_t
        nlh_e_t
        nuh_e_t
        data_size
    end

    methods
        function obj = ZoroDescription()
            % Constructor - initialize the object if needed
        end

        function obj = make_consistent(obj, dims)

            [nw, ~] = size(obj.W_mat);
            obj.nw = nw;
            if isempty(obj.unc_jac_G_mat)
                obj.unc_jac_G_mat = eye(obj.nw);
            end

            % infer np if not set, from Sigma_p_mat if available
            if isempty(obj.np)
                if ~isempty(obj.Sigma_p_mat)
                    obj.np = size(obj.Sigma_p_mat, 1);     % Sigma_p is np x np
                else
                    obj.np = 0;
                end
            end

            obj.nlbx_t = numel(obj.idx_lbx_t);
            obj.nubx_t = numel(obj.idx_ubx_t);
            obj.nlbx_e_t = numel(obj.idx_lbx_e_t);
            obj.nubx_e_t = numel(obj.idx_ubx_e_t);
            obj.nlbu_t = numel(obj.idx_lbu_t);
            obj.nubu_t = numel(obj.idx_ubu_t);
            obj.nlg_t = numel(obj.idx_lg_t);
            obj.nug_t = numel(obj.idx_ug_t);
            obj.nlg_e_t = numel(obj.idx_lg_e_t);
            obj.nug_e_t = numel(obj.idx_ug_e_t);
            obj.nlh_t = numel(obj.idx_lh_t);
            obj.nuh_t = numel(obj.idx_uh_t);
            obj.nlh_e_t = numel(obj.idx_lh_e_t);
            obj.nuh_e_t = numel(obj.idx_uh_e_t);

            % consistency checks for streaming flags
            if obj.input_P0_diag && obj.input_P0
                error('Only one of input_P0_diag and input_P0 can be True. Note: input_P0 defaults to true');
            end
            if obj.input_Sigma_p_diag && obj.input_Sigma_p
                error('Only one of input_Sigma_p_diag and input_Sigma_p can be True');
            end

            FEEDBACK_OPTIMIZATION_MODES = {'CONSTANT_FEEDBACK', 'RICCATI_CONSTANT_COST', 'RICCATI_BARRIER_1', 'RICCATI_BARRIER_2'};

            if ~ismember(obj.feedback_optimization_mode, FEEDBACK_OPTIMIZATION_MODES)
                error('feedback_optimization_mode should be in %s, got %s.', strjoin(FEEDBACK_OPTIMIZATION_MODES, ', '), obj.feedback_optimization_mode);
            end

            if ~strcmp(obj.feedback_optimization_mode, 'CONSTANT_FEEDBACK')
                if isempty(obj.riccati_Q_const) || isempty(obj.riccati_R_const) || isempty(obj.riccati_S_const)
                    error('riccati_Q_const, riccati_R_const, riccati_S_const should not be empty when feedback_optimization_mode ~= CONSTANT_FEEDBACK.');
                end

                if ~isequal(size(obj.riccati_Q_const), [dims.nx, dims.nx])
                    error('The shape of riccati_Q_const should be [nx nx].');
                end
                if ~isequal(size(obj.riccati_R_const), [dims.nu, dims.nu])
                    error('The shape of riccati_R_const should be [nu nu].');
                end
                if ~isequal(size(obj.riccati_S_const), [dims.nu, dims.nx])
                    error('The shape of riccati_S_const should be [nu nx].');
                end

                if isempty(obj.riccati_Q_const_e)
                    obj.riccati_Q_const_e = obj.riccati_Q_const;
                end
                if ~isequal(size(obj.riccati_Q_const_e), [dims.nx, dims.nx])
                    error('The shape of riccati_Q_const_e should be [nx nx].');
                end

            end

            PARAMETER_UNCERTAINTY_MODES = {'CONSTANT', 'IID', 'NONE'};

            if ~ismember(obj.parameter_uncertainty_mode, PARAMETER_UNCERTAINTY_MODES)
                error('parameter_uncertainty_mode should be in %s, got %s.', ...
                    strjoin(PARAMETER_UNCERTAINTY_MODES, ', '), obj.parameter_uncertainty_mode);
            end

            data_size = 0;
            % Print input note:
            fprintf('\nThe data of the generated custom update function consists of the concatenation of:\n');
            i_component = 1;

            % P0 streaming
            if obj.input_P0_diag
                size_i = dims.nx;
                fprintf('%d) input: diag(P0), size: [nx] = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end
            if obj.input_P0
                size_i = dims.nx * dims.nx;
                fprintf('%d) input: P0; full matrix in column-major format, size: [nx*nx] = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end

            % W streaming
            if obj.input_W_diag
                size_i = obj.nw;
                fprintf('%d) input: diag(W), size: [nw] = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end

            % stage-wise additive W_gp^k
            if obj.input_W_add_diag
                size_i = dims.N * obj.nw;
                fprintf('%d) input: concatenation of diag(W_gp^k) for i=0,...,N-1, size: [N * nw] = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end

            % Sigma_p streaming (for S_p term)
            if obj.input_Sigma_p_diag && obj.np > 0
                size_i = obj.np;
                fprintf('%d) input: diag(Sigma_p), size: [np] = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end
            if obj.input_Sigma_p && obj.np > 0
                size_i = obj.np * obj.np;
                fprintf('%d) input: Sigma_p; full matrix in column-major format, size: [np*np] = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end

            % P^k output matrices: output only, not part of payload
            if obj.output_P_matrices
                size_i = dims.nx * dims.nx * (dims.N+1);
                fprintf('%d) output: concatenation of colmaj(P^k) for i=0,...,N, size: [nx*nx*(N+1)] = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end

            if obj.output_riccati_t
                size_i = 1;
                fprintf('%d) output: concatenation of riccati_time = %d\n', i_component, size_i);
                i_component = i_component + 1;
                data_size = data_size + size_i;
            end

            obj.data_size = data_size;
            fprintf('\n');
        end

        function s = to_struct(self)
            if exist('properties')
                publicProperties = eval('properties(self)');
            else
                publicProperties = fieldnames(self);
            end
            s = struct();
            for fi = 1:numel(publicProperties)
                s.(publicProperties{fi}) = self.(publicProperties{fi});
            end
        end

        function s = convert_to_struct_for_json_dump(self, N)
            s = self.to_struct();
            s = prepare_struct_for_json_dump(s, {

                'idx_lbx_t', 'idx_ubx_t', 'idx_lbx_e_t', 'idx_ubx_e_t', ...
                'idx_lbu_t', 'idx_ubu_t', 'idx_lg_t', 'idx_ug_t', ...
                'idx_lg_e_t', 'idx_ug_e_t', 'idx_lh_t', 'idx_uh_t', ...
                'idx_lh_e_t', 'idx_uh_e_t'}, {
                    'fdbk_K_mat', 'unc_jac_G_mat', 'P0_mat', 'W_mat', 'Sigma_p_mat'});

        end
    end
end

classdef ZoroDescription < handle
    properties
        backoff_scaling_gamma = 1.0
        fdbk_K_mat = []
        unc_jac_G_mat = []
        P0_mat = []
        W_mat = []

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

        input_P0_diag = false
        input_P0 = true
        input_W_diag = false
        input_W_add_diag = false

        output_P_matrices = false

    % properties (Access = private)
    % kind of private, but need to be dumped to json
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
    end

    methods
        function obj = ZoroDescription()
            % Constructor - initialize the object if needed
        end

        function obj = process(obj)
            [nw, ~] = size(obj.W_mat);
            obj.nw = nw;
            if isempty(obj.unc_jac_G_mat)
                obj.unc_jac_G_mat = eye(obj.nw);
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

            if obj.input_P0_diag && obj.input_P0
                error('Only one of input_P0_diag and input_P0 can be True');
            end

            % Print input note:
            fprintf('\nThe data of the generated custom update function consists of the concatenation of:\n');
            i_component = 1;
            if obj.input_P0_diag
                fprintf('%d) input: diag(P0)\n', i_component);
                i_component = i_component + 1;
            end
            if obj.input_P0
                fprintf('%d) input: P0; full matrix in column-major format\n', i_component);
                i_component = i_component + 1;
            end
            if obj.input_W_diag
                fprintf('%d) input: diag(W)\n', i_component);
                i_component = i_component + 1;
            end
            if obj.input_W_add_diag
                fprintf('%d) input: concatenation of diag(W_gp^k) for i=0,...,N-1\n', i_component);
                i_component = i_component + 1;
            end
            if obj.output_P_matrices
                fprintf('%d) output: concatenation of colmaj(P^k) for i=0,...,N\n', i_component);
            end
            fprintf('\n');
        end

        function s = struct(self)
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
            s = self.struct();
            s = prepare_struct_for_json_dump(s, {
                'idx_lbx_t', 'idx_ubx_t', 'idx_lbx_e_t', 'idx_ubx_e_t', 'idx_lbu_t', 'idx_ubu_t', 'idx_lg_t', 'idx_ug_t', 'idx_lg_e_t', 'idx_ug_e_t', 'idx_lh_t', 'idx_uh_t', 'idx_lh_e_t', 'idx_uh_e_t'}, {
                    'fdbk_K_mat', 'unc_jac_G_mat', 'P0_mat', 'W_mat'});
        end
    end
end

%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

classdef AcadosOcp < handle
    properties
        dims
        cost
        constraints
        solver_options
        model
        parameter_values % initial value of the parameter
        acados_include_path
        acados_lib_path
        problem_class
        simulink_opts
        cython_include_dirs
        code_export_directory
        json_file
        shared_lib_ext
        name
    end
    methods
        function obj = AcadosOcp()
            obj.dims = AcadosOcpDims();
            obj.cost = AcadosOcpCost();
            obj.constraints = AcadosOcpConstraints();
            obj.solver_options = AcadosOcpOptions();
            obj.model = AcadosModel();
            obj.acados_include_path = [];
            obj.acados_lib_path = [];
            obj.parameter_values = [];
            obj.problem_class = 'OCP';
            obj.simulink_opts = [];
            obj.cython_include_dirs = [];
            obj.json_file = 'acados_ocp_nlp.json';
            obj.shared_lib_ext = '.so';
            obj.name = 'ocp';
            if ismac()
                obj.shared_lib_ext = '.dylib';
            end
            obj.code_export_directory = 'c_generated_code';
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

        function make_consistent(self)
            self.model.make_consistent(self.dims);

            N = self.dims.N;

            model = self.model;
            dims = self.dims;
            cost = self.cost;
            constraints = self.constraints;
            opts = self.solver_options;

            self.detect_cost_and_constraints();
            self.name = model.name;

            % parameters
            if isempty(self.parameter_values)
                warning(['self.parameter_values are not set.', ...
                            10 'Using zeros(np,1) by default.' 10 'You can update them later using set().']);
                self.parameter_values = zeros(self.dims.np,1);
            elseif length(self.parameter_values) ~= self.dims.np
                error(['parameters_values has the wrong shape. Expected: ' num2str(self.dims.np)])
            end

            %% cost
            % initial
            if strcmp(cost.cost_type_0, 'LINEAR_LS')
                if ~isempty(cost.W_0) && ~isempty(cost.Vx_0) && ~isempty(cost.Vu_0)
                    ny = length(cost.W_0);

                    if isempty(cost.yref_0)
                        warning(['yref_0 not defined provided.' 10 'Using zeros(ny_0,1) by default.']);
                        self.cost.yref_0 = zeros(ny,1);
                    end
                    if ny ~= size(cost.Vx_0, 1) || ny ~= size(cost.Vu_0, 1) || ny ~= size(cost.yref_0, 1)
                        error('inconsistent dimension ny, regarding W_0, Vx_0, Vu_0, yref_0.');
                    end
                else
                    error('setting linear least square cost: need W_0, Vx_0, Vu_0, at least one missing.')
                end
                dims.ny_0 = ny;
            elseif strcmp(cost.cost_type_0, 'NONLINEAR_LS')
                if ~isempty(cost.W_0) && ~isempty(model.cost_y_expr_0)
                    ny = length(cost.W_0);
                    if isempty(cost.yref_0)
                        warning(['yref_0 not defined provided.' 10 'Using zeros(ny_0,1) by default.']);
                        self.cost.yref_0 = zeros(ny,1);
                    end
                    if ny ~= length(model.cost_y_expr_0) || ny ~= size(cost.yref_0, 1)
                        error('inconsistent dimension ny_0, regarding W_0, cost_y_expr_0, yref_0.');
                    end
                else
                    error('setting nonlinear least square cost: need W_0, cost_y_expr_0, at least one missing.')
                end
                dims.ny_0 = ny;
            end

            % path
            if strcmp(cost.cost_type, 'LINEAR_LS')
                if ~isempty(cost.W) && ~isempty(cost.Vx) && ~isempty(cost.Vu)
                    ny = length(cost.W);
                    if isempty(cost.yref)
                        warning(['yref not defined provided.' 10 'Using zeros(ny,1) by default.']);
                        self.cost.yref = zeros(ny,1);
                    end
                    if ny ~= size(cost.Vx, 1) || ny ~= size(cost.Vu, 1) || ny ~= size(cost.yref, 1)
                        error('inconsistent dimension ny, regarding W, Vx, Vu, yref.');
                    end
                else
                    error('setting linear least square cost: need W, Vx, Vu, at least one missing.')
                end
                dims.ny = ny;
            elseif strcmp(cost.cost_type, 'NONLINEAR_LS')
                if ~isempty(cost.W) && ~isempty(model.cost_y_expr)
                    ny = length(cost.W);
                    if isempty(cost.yref)
                        warning(['yref not defined provided.' 10 'Using zeros(ny,1) by default.']);
                        self.cost.yref = zeros(ny,1);
                    end
                    if ny ~= length(model.cost_y_expr) || ny ~= size(cost.yref, 1)
                        error('inconsistent dimension ny, regarding W, cost_y_expr, yref.');
                    end
                else
                    error('setting nonlinear least square cost: need W, cost_y_expr, at least one missing.')
                end
                dims.ny = ny;
            end

            % terminal
            if strcmp(cost.cost_type_e, 'LINEAR_LS')
                if ~isempty(cost.W_e) && ~isempty(cost.Vx_e)
                    ny_e = length(cost.W_e);
                    if isempty(cost.yref_e)
                        warning(['yref_e not defined provided.' 10 'Using zeros(ny_e,1) by default.']);
                        self.cost.yref_e = zeros(ny_e,1);
                    end
                    if ny_e ~= size(cost.Vx_e, 1) || ny_e ~= size(cost.yref_e, 1)
                        error('inconsistent dimension ny_e, regarding W_e, Vx_e, yref_e');
                    end
                elseif ~~isempty(cost.W_e) && ~~isempty(cost.Vx_e)
                    ny_e = 0;
                    warning('Fields W_e and Vx_e not provided. Using empty ls terminal cost.')
                else
                    error('setting linear least square cost: need W_e, Vx_e, at least one missing.')
                end
                dims.ny_e = ny_e;
            elseif strcmp(cost.cost_type_e, 'NONLINEAR_LS')
                if ~isempty(cost.W_e) && ~isempty(model.cost_y_expr_e)
                    ny_e = length(cost.W_e);
                    if isempty(cost.yref_e)
                        warning(['yref_e not defined provided.' 10 'Using zeros(ny_e,1) by default.']);
                        self.cost.yref_e = zeros(ny_e,1);
                    end
                    if ny_e ~= length(model.cost_y_expr_e) || ny_e ~= size(cost.yref_e, 1)
                        error('inconsistent dimension ny_e, regarding W_e, cost_y_expr_e, yref_e.');
                    end
                else
                    error('setting nonlinear least square cost: need W_e, cost_y_expr_e, at least one missing.')
                end
                dims.ny_e = ny_e;
            end


            %% constraints
            % initial
            if ~isempty(constraints.x0)

                if length(constraints.x0) ~= dims.nx
                    error('inconsistent constraint x0, regarding nx.');
                end
                if isempty(constraints.idxbx_0) && isempty(constraints.lbx_0) && isempty(constraints.ubx_0) && isempty(constraints.idxbxe_0)
                    constraints.idxbx_0 = 0:dims.nx-1;
                    constraints.idxbxe_0 = 0:dims.nx-1;
                    constraints.lbx_0 = constraints.x0;
                    constraints.ubx_0 = constraints.x0;
                    constraints.has_x0 = true;
                else
                    error('If constraint x0 is defined, bounds on x0 should not be defined, check properties lbx_0, ubx_0, idxbx_0, idxbxe_0.');
                end
            end
            if ~isempty(constraints.idxbx_0) && ~isempty(constraints.lbx_0) && ~isempty(constraints.ubx_0)
                nbx_0 = length(constraints.lbx_0);
                if nbx_0 ~= length(constraints.ubx_0) || nbx_0 ~= length(constraints.idxbx_0)
                    error('inconsistent dimension nbx_0, regarding Jbx_0, lbx_0, ubx_0.');
                end
            elseif ~isempty(constraints.idxbx_0) || ~isempty(constraints.lbx_0) || ~isempty(constraints.ubx_0)
                error('setting bounds on x: need Jbx_0, lbx_0, ubx_0, at least one missing.');
            else
                % no initial state constraint
                disp("OCP without constraints on initial state detected.")
                nbx_0 = 0;
            end
            dims.nbx_0 = nbx_0;

            if ~isempty(constraints.idxbxe_0)
                dims.nbxe_0 = length(constraints.idxbxe_0);
            else
                % no equalities on initial state.
                constraints.idxbxe_0 = [];
                dims.nbxe_0 = 0;
            end

            % path
            if ~isempty(constraints.idxbx) && ~isempty(constraints.lbx) && ~isempty(constraints.ubx)
                nbx = length(constraints.lbx);
                if nbx ~= length(constraints.ubx) || nbx ~= length(constraints.idxbx)
                    error('inconsistent dimension nbx, regarding Jbx, lbx, ubx.');
                end
            elseif ~isempty(constraints.idxbx) || ~isempty(constraints.lbx) || ~isempty(constraints.ubx)
                error('setting bounds on x: need Jbx, lbx, ubx, at least one missing.');
            else
                nbx = 0;
            end
            dims.nbx = nbx;

            if ~isempty(constraints.idxbu) && ~isempty(constraints.lbu) && ~isempty(constraints.ubu)
                nbu = length(constraints.lbu);
                if nbu ~= length(constraints.ubu) || nbu ~= length(constraints.idxbu)
                    error('inconsistent dimension nbu, regarding Jbu, lbu, ubu.');
                end
            elseif ~isempty(constraints.idxbu) || ~isempty(constraints.lbu) || ~isempty(constraints.ubu)
                error('setting bounds on u: need Jbu, lbu, ubu, at least one missing.');
            else
                nbu = 0;
            end
            dims.nbu = nbu;

            if ~isempty(constraints.C) && ~isempty(constraints.D) && ...
               ~isempty(constraints.lg) && ~isempty(constraints.ug)
                ng = length(constraints.lg);
                if ng ~= length(constraints.ug) || ng ~= size(constraints.C, 1) || ng ~= size(constraints.D, 1)
                    error('inconsistent dimension ng, regarding C, D, lg, ug.');
                end
            elseif ~isempty(constraints.C) || ~isempty(constraints.D) || ...
                   ~isempty(constraints.lg) || ~isempty(constraints.ug)
                error('setting general linear constraints: need C, D, lg, ug, at least one missing.');
            else
                ng = 0;
            end
            dims.ng = ng;

            if ~isempty(model.con_h_expr) && ...
                     ~isempty(constraints.lh) && ~isempty(constraints.uh)
                nh = length(constraints.lh);
                if nh ~= length(constraints.uh) || nh ~= length(model.con_h_expr)
                    error('inconsistent dimension nh, regarding expr_h, lh, uh.');
                end
            elseif ~isempty(model.con_h_expr) || ...
                   ~isempty(constraints.lh) || ~isempty(constraints.uh)
                error('setting external constraint function h: need expr_h, lh, uh at least one missing.');
            else
                nh = 0;
            end
            dims.nh = nh;

            if ~isempty(model.con_h_expr_0) && ...
                    ~isempty(constraints.lh_0) && ~isempty(constraints.uh_0)
            nh_0 = length(constraints.lh_0);
            if nh_0 ~= length(constraints.uh_0) || nh_0 ~= length(model.con_h_expr_0)
                error('inconsistent dimension nh_0, regarding expr_h_0, lh_0, uh_0.');
            end
            elseif ~isempty(model.con_h_expr_0) || ...
                ~isempty(constraints.lh_0) || ~isempty(constraints.uh_0)
            error('setting external constraint function h: need expr_h_0, lh_0, uh_0 at least one missing.');
            else
                nh_0 = 0;
            end
            dims.nh_0 = nh_0;

            % terminal
            if ~isempty(constraints.idxbx_e) && ~isempty(constraints.lbx_e) && ~isempty(constraints.ubx_e)
                nbx_e = length(constraints.lbx_e);
                if nbx_e ~= length(constraints.ubx_e) || nbx_e ~= length(constraints.idxbx_e)
                    error('inconsistent dimension nbx_e, regarding Jbx_e, lbx_e, ubx_e.');
                end
            elseif ~isempty(constraints.idxbx_e) || ~isempty(constraints.lbx_e) || ~isempty(constraints.ubx_e)
                error('setting bounds on x: need Jbx_e, lbx_e, ubx_e, at least one missing.');
            else
                nbx_e = 0;
            end
            dims.nbx_e = nbx_e;

            if ~isempty(constraints.C_e) && ...
               ~isempty(constraints.lg_e) && ~isempty(constraints.ug_e)
                ng_e = length(constraints.lg_e);
                if ng_e ~= length(constraints.ug_e) || ng_e ~= size(constraints.C_e, 1)
                    error('inconsistent dimension ng_e, regarding C_e, lg_e, ug_e.');
                end
            elseif ~isempty(constraints.C_e) || ...
                   ~isempty(constraints.lg_e) || ~isempty(constraints.ug_e)
                error('setting general linear constraints: need C_e, lg_e, ug_e, at least one missing.');
            else
                ng_e = 0;
            end
            dims.ng_e = ng_e;

            if ~isempty(model.con_h_expr_e) && ...
                     ~isempty(constraints.lh_e) && ~isempty(constraints.uh_e)
                nh_e = length(constraints.lh_e);
                if nh_e ~= length(constraints.uh_e) || nh_e ~= length(model.con_h_expr_e)
                    error('inconsistent dimension nh_e, regarding expr_h_e, lh_e, uh_e.');
                end
            elseif ~isempty(model.con_h_expr_e) || ...
                   ~isempty(constraints.lh_e) || ~isempty(constraints.uh_e)
                error('setting external constraint function h: need expr_h_e, lh_e, uh_e at least one missing.');
            else
                nh_e = 0;
            end
            dims.nh_e = nh_e;

            %% slack dimensions
            if ~isempty(constraints.idxsbx)
                nsbx = length(constraints.idxsbx);
            else
                nsbx = 0;
            end

            if ~isempty(constraints.idxsbu)
                nsbu = length(constraints.idxsbu);
            else
                nsbu = 0;
            end

            if ~isempty(constraints.idxsg)
                nsg = length(constraints.idxsg);
            else
                nsg = 0;
            end
            if ~isempty(constraints.idxsh)
                nsh = length(constraints.idxsh);
            else
                nsh = 0;
            end
            if ~isempty(constraints.idxsphi)
                nsphi = length(constraints.idxsphi);
            else
                nsphi = 0;
            end

            ns = nsbx + nsbu + nsg + nsh + nsphi;
            wrong_field = '';
            if ~isempty(cost.Zl) && ~all(size(cost.Zl) == [ns, 1])
                wrong_field = 'Zl';
                dim = size(cost.Zl);
            elseif ~isempty(cost.Zu) && ~all(size(cost.Zu) == [ns, 1])
                wrong_field = 'Zu';
                dim = size(cost.Zu);
            elseif ~isempty(cost.zl) && ~all(size(cost.zl) == [ns, 1])
                wrong_field = 'zl';
                dim = size(cost.zl);
            elseif ~isempty(cost.zu) && ~all(size(cost.zu) == [ns, 1])
                wrong_field = 'zu';
                dim = size(cost.zu);
            end

            if ~strcmp(wrong_field, '')
                error(['Inconsistent size for field ', wrong_field, ' with dimension ', num2str(dim),...
                      '. Detected ns = ', num2str(ns), ' = nsbx + nsbu + nsg + nsh + nsphi.',...
                      ' With nsbx = ', num2str(nsbx), ', nsbu = ', num2str(nsbu), ' nsg = ', num2str(nsg),...
                      ' nsh = ', num2str(nsh), ', nsphi = ', num2str(nsphi), '.'])
            end

            constraints.lsbu = zeros(nsbu, 1);
            constraints.usbu = zeros(nsbu, 1);
            constraints.lsbx = zeros(nsbx, 1);
            constraints.usbx = zeros(nsbx, 1);
            constraints.lsh = zeros(nsh, 1);
            constraints.ush = zeros(nsh, 1);
            constraints.lsphi = zeros(nsphi, 1);
            constraints.usphi = zeros(nsphi, 1);

            dims.ns = ns;
            dims.nsbx = nsbx;
            dims.nsbu = nsbu;
            dims.nsg = nsg;
            dims.nsh = nsh;
            dims.nsphi = nsphi;

            % slacks at initial stage
            if ~isempty(constraints.idxsh_0)
                nsh_0 = length(constraints.idxsh_0);
            else
                nsh_0 = 0;
            end
            if ~isempty(constraints.idxsphi_0)
                nsphi_0 = length(constraints.idxsphi_0);
            else
                nsphi_0 = 0;
            end

            ns_0 = nsbu + nsg + nsh_0 + nsphi_0;
            wrong_field = '';
            if ~isempty(cost.Zl_0) && ~all(size(cost.Zl_0) == [ns_0, 1])
                wrong_field = 'Zl_0';
                dim = size(cost.Zl_0);
            elseif ~isempty(cost.Zu_0) && ~all(size(cost.Zu_0) == [ns_0, 1])
                wrong_field = 'Zu_0';
                dim = size(cost.Zu_0);
            elseif ~isempty(cost.zl_0) && ~all(size(cost.zl_0) == [ns_0, 1])
                wrong_field = 'zl_0';
                dim = size(cost.zl_0);
            elseif ~isempty(cost.zu_0) && ~all(size(cost.zu_0) == [ns_0, 1])
                wrong_field = 'zu_0';
                dim = size(cost.zu_0);
            end

            if ~strcmp(wrong_field, '')
                error(['Inconsistent size for field', wrong_field, ' with dimension ', num2str(dim),...
                        '. Detected ns_0 = ', num2str(ns_0), ' = nsbu + nsg + nsh_0 + nsphi_0.',...
                        ' With nsg = ', num2str(nsg), ' nsh_0 = ', num2str(nsh_0), ', nsphi_0 = ', num2str(nsphi_0), '.'])
            end

            constraints.lsh_0 = zeros(nsh_0, 1);
            constraints.ush_0 = zeros(nsh_0, 1);
            constraints.lsphi_0 = zeros(nsphi_0, 1);
            constraints.usphi_0 = zeros(nsphi_0, 1);

            dims.ns_0 = ns_0;
            dims.nsh_0 = nsh_0;
            dims.nsphi_0 = nsphi_0;

            %% terminal slack dimensions
            if ~isempty(constraints.idxsbx_e)
                nsbx_e = length(constraints.idxsbx_e);
            else
                nsbx_e = 0;
            end

            if ~isempty(constraints.idxsg_e)
                nsg_e = length(constraints.idxsg_e);
            else
                nsg_e = 0;
            end
            if ~isempty(constraints.idxsh_e)
                nsh_e = length(constraints.idxsh_e);
            else
                nsh_e = 0;
            end
            if ~isempty(constraints.idxsphi_e)
                nsphi_e = length(constraints.idxsphi_e);
            else
                nsphi_e = 0;
            end

            ns_e = nsbx_e + nsg_e + nsh_e + nsphi_e;
            wrong_field = '';
            if ~isempty(cost.Zl_e) && ~all(size(cost.Zl_e) == [ns_e, 1])
                wrong_field = 'Zl_e';
                dim = size(cost.Zl_e);
            elseif ~isempty(cost.Zu_e) && ~all(size(cost.Zu_e) == [ns_e, 1])
                wrong_field = 'Zu_e';
                dim = size(cost.Zu_e);
            elseif ~isempty(cost.zl_e) && ~all(size(cost.zl_e) == [ns_e, 1])
                wrong_field = 'zl_e';
                dim = size(cost.zl_e);
            elseif ~isempty(cost.zu_e) && ~all(size(cost.zu_e) == [ns_e, 1])
                wrong_field = 'zu_e';
                dim = size(cost.zu_e);
            end

            if ~strcmp(wrong_field, '')
                error(['Inconsistent size for field', wrong_field, ' with dimension ', num2str(dim),...
                        '. Detected ns_e = ', num2str(ns_e), ' = nsbx_e + nsg_e + nsh_e + nsphi_e.',...
                        ' With nsbx_e = ', num2str(nsbx_e), ' nsg_e = ', num2str(nsg_e),...
                        ' nsh_e = ', num2str(nsh_e), ', nsphi_e = ', num2str(nsphi_e), '.'])
            end


            constraints.lsbx_e = zeros(nsbx_e, 1);
            constraints.usbx_e = zeros(nsbx_e, 1);
            constraints.lsh_e = zeros(nsh_e, 1);
            constraints.ush_e = zeros(nsh_e, 1);
            constraints.lsphi_e = zeros(nsphi_e, 1);
            constraints.usphi_e = zeros(nsphi_e, 1);

            dims.ns_e = ns_e;
            dims.nsbx_e = nsbx_e;
            dims.nsg_e = nsg_e;
            dims.nsh_e = nsh_e;
            dims.nsphi_e = nsphi_e;

            % shooting nodes -> time_steps
            N = self.dims.N;
            if length(opts.tf) ~= 1 || opts.tf < 0
                error('time horizon tf should be a nonnegative number');
            end

            if ~isempty(opts.shooting_nodes)
                if N + 1 ~= length(opts.shooting_nodes)
                    error('inconsistent dimension N regarding shooting nodes.');
                end
                for i=1:N
                    opts.time_steps(i) = opts.shooting_nodes(i+1) - opts.shooting_nodes(i);
                end
                sum_time_steps = sum(opts.time_steps);
                if abs((sum_time_steps - opts.tf) / opts.tf) > 1e-14
                    warning('shooting nodes are not consistent with time horizon tf, rescaling automatically');
                    opts.time_steps = opts.time_steps * opts.tf / sum_time_steps;
                end
            elseif ~isempty(opts.time_steps)
                if N ~= length(opts.time_steps)
                    error('inconsistent dimension N regarding time steps.');
                end
                sum_time_steps = sum(opts.time_steps);
                if abs((sum_time_steps - opts.tf) / opts.tf) > 1e-14
                    error(['time steps are not consistent with time horizon tf, ', ...
                        'got tf = ' num2str(opts.tf) '; sum(time_steps) = ' num2str(sum_time_steps) '.']);
                end
                % just to have them available, e.g. for plotting;
                opts.shooting_nodes = zeros(N+1, 1);
                for i = 1:N
                    opts.shooting_nodes(i+1) = sum(opts.time_steps(1:i));
                end
            else
                opts.time_steps = opts.tf/N * ones(N,1);
            end
            if any(opts.time_steps < 0)
                error(['ocp discretization: time_steps between shooting nodes must all be > 0', ...
                    ' got: ' num2str(opts.time_steps)])
            end
            if ~isempty(opts.sim_method_num_stages)
                if(strcmp(opts.integrator_type, "ERK"))
                    if (any(opts.sim_method_num_stages < 1) || any(opts.sim_method_num_stages > 4))
                        error(['ERK: num_stages = ', num2str(opts.sim_method_num_stages) ' not available. Only number of stages = {1,2,3,4} implemented!']);
                    end
                end
            end

            % set integrator time automatically
            opts.Tsim = opts.time_steps(1);

            % qpdunes
            if ~isempty(strfind(opts.qp_solver,'qpdunes'))
                constraints.idxbxe_0 = [];
                dims.nbxe_0 = 0;
            end

            % options sanity checks
            if length(self.solver_options.sim_method_num_steps) ~= N
                self.solver_options.sim_method_num_steps = self.solver_options.sim_method_num_steps * ones(1, N);
            end
            if length(self.solver_options.sim_method_num_stages) ~= N
                self.solver_options.sim_method_num_stages = self.solver_options.sim_method_num_stages * ones(1, N);
            end
            if length(self.solver_options.sim_method_jac_reuse) ~= N
                self.solver_options.sim_method_jac_reuse = self.solver_options.sim_method_jac_reuse * ones(1, N);
            end

            if strcmp(self.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_OSQP") || ...
                strcmp(self.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_HPMPC") || ...
                strcmp(self.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_QPDUNES") || ...
                strcmp(self.solver_options.nlp_solver_type, "PARTIAL_CONDENSING_OOQP")
                if self.dims.ns > 0 || self.dims.ns_e > 0
                    error(['selected QP solver ', self.solver_options.nlp_solver_type, ' does not support soft constraints (yet).'])
                end
            end
        end

        function [] = detect_cost_and_constraints(self)
            % detect cost type
            stage_types = {'initial', 'path', 'terminal'};
            cost_types = {self.cost.cost_type_0, self.cost.cost_type, self.cost.cost_type_e};

            for n=1:3
                if strcmp(cost_types{n}, 'AUTO')
                    detect_cost_type(self.model, self.cost, self.dims, stage_types{n});
                end
            end

            % if initial is empty, copy path cost
            % TODO: move this to make_consistent? should happen way before?
            if isempty(cost_types{1})
                warning("cost_type_0 not set, using path cost");
                self.cost.cost_type_0 = self.cost.cost_type;
                if (strcmp(self.cost.cost_type, 'LINEAR_LS'))
                    self.cost.Vx_0 = self.cost.Vx;
                    self.cost.Vu_0 = self.cost.Vu;
                    self.cost.Vz_0 = self.cost.Vz;
                elseif (strcmp(self.cost.cost_type, 'NONLINEAR_LS'))
                    self.model.cost_y_expr_0 = self.model.cost_y_expr;
                elseif (strcmp(self.cost.cost_type, 'EXTERNAL'))
                    self.cost.cost_ext_fun_type_0 = self.cost.cost_ext_fun_type;
                    if strcmp(self.cost.cost_ext_fun_type_0, 'casadi')
                        self.model.cost_expr_ext_cost_0 = self.model.cost_expr_ext_cost;
                        self.model.cost_expr_ext_cost_custom_hess_0 = self.model.cost_expr_ext_cost_custom_hess;
                    else % generic
                        self.cost.cost_source_ext_cost_0 = self.cost.cost_source_ext_cost;
                        self.cost.cost_function_ext_cost_0 = self.cost.cost_function_ext_cost;
                    end
                end
                if (strcmp(self.cost.cost_type, 'LINEAR_LS')) || (strcmp(self.cost.cost_type, 'NONLINEAR_LS'))
                    self.cost.W_0 = self.cost.W;
                    self.cost.yref_0 = self.cost.yref;
                    self.dims.ny_0 = self.dims.ny;
                end
            end

            % detect constraint structure
            constraint_types = {self.constraints.constr_type_0, self.constraints.constr_type, self.constraints.constr_type_e};
            for n=1:3
                if strcmp(constraint_types{n}, 'AUTO')
                    detect_constr(self.model, self.constraints, stage_types{n});
                end
            end
        end
    end
end


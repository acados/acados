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
        p_global_values % initial value of the parameter
        acados_include_path
        acados_lib_path
        problem_class
        simulink_opts
        cython_include_dirs
        code_export_directory
        json_file
        shared_lib_ext
        name
        zoro_description
        external_function_files_ocp
        external_function_files_model
    end
    methods
        function obj = AcadosOcp()
            obj.dims = AcadosOcpDims();
            obj.cost = AcadosOcpCost();
            obj.constraints = AcadosOcpConstraints();
            obj.solver_options = AcadosOcpOptions();
            obj.model = AcadosModel();

            obj.parameter_values = [];
            obj.p_global_values = [];
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

            % set include and lib path
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            obj.acados_include_path = [acados_folder, '/include'];
            obj.acados_lib_path = [acados_folder, '/lib'];
            obj.zoro_description = [];
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

        function make_consistent_cost_initial(self)
            cost = self.cost;
            dims = self.dims;
            model = self.model;
            if self.solver_options.N_horizon == 0
                return
            end
            if strcmp(cost.cost_type_0, 'LINEAR_LS')
                if ~isempty(cost.W_0) && ~isempty(cost.Vx_0) && ~isempty(cost.Vu_0)
                    ny = length(cost.W_0);

                    if isempty(cost.yref_0)
                        warning(['yref_0 not defined provided.' 10 'Using zeros(ny_0,1) by default.']);
                        self.cost.yref_0 = zeros(ny,1);
                    end
                    if ny ~= size(cost.Vx_0, 1) || ny ~= size(cost.Vu_0, 1) || ny ~= size(cost.yref_0, 1)
                        error('inconsistent dimension ny_0, regarding W_0, Vx_0, Vu_0, yref_0.');
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
        end

        function make_consistent_cost_path(self)
            cost = self.cost;
            dims = self.dims;
            model = self.model;
            if self.solver_options.N_horizon == 0
                return
            end
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
        end

        function make_consistent_cost_terminal(self)
            cost = self.cost;
            dims = self.dims;
            model = self.model;

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
        end

        function make_consistent_constraints_initial(self)
            dims = self.dims;
            constraints = self.constraints;
            model = self.model;
            if self.solver_options.N_horizon == 0
                dims.nbxe_0 = 0;
                return
            end

            if ~isempty(constraints.idxbx_0) && ~isempty(constraints.lbx_0) && ~isempty(constraints.ubx_0)
                nbx_0 = length(constraints.lbx_0);
                if nbx_0 ~= length(constraints.ubx_0) || nbx_0 ~= length(constraints.idxbx_0)
                    error('inconsistent dimension nbx_0, regarding idxbx_0, lbx_0, ubx_0.');
                end
                if min(constraints.idxbx_0) < 0 || max(constraints.idxbx_0) > (dims.nx-1)
                    error(['idxbx_0 should contain (zero-based) indices between 0 and ', num2str(dims.nx-1)])
                end
            elseif ~isempty(constraints.idxbx_0) || ~isempty(constraints.lbx_0) || ~isempty(constraints.ubx_0)
                error('setting bounds on x: need idxbx_0, lbx_0, ubx_0, at least one missing.');
            else
                % no initial state constraint
                disp("OCP without constraints on initial state detected.")
                nbx_0 = 0;
            end
            dims.nbx_0 = nbx_0;

            dims.nbxe_0 = length(constraints.idxbxe_0);

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
        end

        function make_consistent_constraints_path(self)
            dims = self.dims;
            constraints = self.constraints;
            model = self.model;
            if self.solver_options.N_horizon == 0
                return
            end

            if ~isempty(constraints.idxbx) && ~isempty(constraints.lbx) && ~isempty(constraints.ubx)
                nbx = length(constraints.lbx);
                if nbx ~= length(constraints.ubx) || nbx ~= length(constraints.idxbx)
                    error('inconsistent dimension nbx, regarding idxbx, lbx, ubx.');
                end
                if min(constraints.idxbx) < 0 || max(constraints.idxbx) > (dims.nx-1)
                    error(['idxbx should contain (zero-based) indices between 0 and ', num2str(dims.nx-1)])
                end
            elseif ~isempty(constraints.idxbx) || ~isempty(constraints.lbx) || ~isempty(constraints.ubx)
                error('setting bounds on x: need idxbx, lbx, ubx, at least one missing.');
            else
                nbx = 0;
            end
            dims.nbx = nbx;

            if ~isempty(constraints.idxbu) && ~isempty(constraints.lbu) && ~isempty(constraints.ubu)
                nbu = length(constraints.lbu);
                if nbu ~= length(constraints.ubu) || nbu ~= length(constraints.idxbu)
                    error('inconsistent dimension nbu, regarding idxbu, lbu, ubu.');
                end
                if min(constraints.idxbu) < 0 || max(constraints.idxbu) > (dims.nu-1)
                    error(['idxbu should contain (zero-based) indices between 0 and ', num2str(dims.nu-1)])
                end
            elseif ~isempty(constraints.idxbu) || ~isempty(constraints.lbu) || ~isempty(constraints.ubu)
                error('setting bounds on u: need idxbu, lbu, ubu, at least one missing.');
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
        end

        function make_consistent_constraints_terminal(self)
            dims = self.dims;
            constraints = self.constraints;
            model = self.model;

            if ~isempty(constraints.idxbx_e) && ~isempty(constraints.lbx_e) && ~isempty(constraints.ubx_e)
                nbx_e = length(constraints.lbx_e);
                if nbx_e ~= length(constraints.ubx_e) || nbx_e ~= length(constraints.idxbx_e)
                    error('inconsistent dimension nbx_e, regarding Jbx_e, lbx_e, ubx_e.');
                end
                if min(constraints.idxbx_e) < 0 || max(constraints.idxbx_e) > (dims.nx-1)
                    error(['idxbx_e should contain (zero-based) indices between 0 and ', num2str(dims.nx-1)])
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
        end

        function make_consistent_slacks_path(self)
            constraints = self.constraints;
            dims = self.dims;
            cost = self.cost;
            if self.solver_options.N_horizon == 0
                return
            end

            nbx = dims.nbx;
            nsbx = length(constraints.idxsbx);
            if nsbx > nbx
                error(['inconsistent dimension nsbx = ', num2str(nsbx), '. Is greater than nbx = ', num2str(nbx), '.']);
            end
            if any(constraints.idxsbx >= nbx)
                error(['idxsbx = [', num2str(constraints.idxsbx(:).'), '] contains value >= nbx = ', num2str(nbx), '.']);
            end
            if any(constraints.idxsbx < 0)
                error(['idxsbx = [', num2str(constraints.idxsbx(:).'), '] contains value < 0.']);
            end

            nsbu = length(constraints.idxsbu);
            nbu = dims.nbu;
            if nsbu > nbu
                error(['inconsistent dimension nsbu = ', num2str(nsbu), '. Is greater than nbu = ', num2str(nbu), '.']);
            end
            if any(constraints.idxsbu >= nbu)
                error(['idxsbu = [', num2str(constraints.idxsbu(:).'), '] contains value >= nbu = ', num2str(nbu), '.']);
            end
            if any(constraints.idxsbu < 0)
                error(['idxsbu = [', num2str(constraints.idxsbu(:).'), '] contains value < 0.']);
            end

            nsg = length(constraints.idxsg);
            ng = dims.ng;
            if nsg > ng
                error(['inconsistent dimension nsg = ', num2str(nsg), '. Is greater than ng = ', num2str(ng), '.']);
            end
            if any(constraints.idxsg >= ng)
                error(['idxsg = [', num2str(constraints.idxsg(:).'), '] contains value >= ng = ', num2str(ng), '.']);
            end
            if any(constraints.idxsg < 0)
                error(['idxsg = [', num2str(constraints.idxsg(:).'), '] contains value < 0.']);
            end

            nsh = length(constraints.idxsh);
            nh = dims.nh;
            if nsh > nh
                error(['inconsistent dimension nsh = ', num2str(nsh), '. Is greater than nh = ', num2str(nh), '.']);
            end
            if any(constraints.idxsh >= nh)
                error(['idxsh = [', num2str(constraints.idxsh(:).'), '] contains value >= nh = ', num2str(nh), '.']);
            end
            if any(constraints.idxsh < 0)
                error(['idxsh = [', num2str(constraints.idxsh(:).'), '] contains value < 0.']);
            end

            nsphi = length(constraints.idxsphi);
            nphi = dims.nphi;
            if nsphi > nphi
                error(['inconsistent dimension nsphi = ', num2str(nsphi), '. Is greater than nphi = ', num2str(nphi), '.']);
            end
            if any(constraints.idxsphi >= nphi)
                error(['idxsphi = [', num2str(constraints.idxsphi(:).'), '] contains value >= nphi = ', num2str(nphi), '.']);
            end
            if any(constraints.idxsphi < 0)
                error(['idxsphi = [', num2str(constraints.idxsphi(:).'), '] contains value < 0.']);
            end

            ns = nsbx + nsbu + nsg + nsh + nsphi;
            wrong_field = '';

            if ns == 0
                expected_shape = [0, 0];
            else
                expected_shape = [ns, 1];
            end

            if ~all(size(cost.Zl) == expected_shape)
                wrong_field = 'Zl';
                dim = size(cost.Zl);
            elseif ~all(size(cost.Zu) == expected_shape)
                wrong_field = 'Zu';
                dim = size(cost.Zu);
            elseif ~all(size(cost.zl) == expected_shape)
                wrong_field = 'zl';
                dim = size(cost.zl);
            elseif ~all(size(cost.zu) == expected_shape)
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
        end

        function make_consistent_slacks_initial(self)
            constraints = self.constraints;
            dims = self.dims;
            cost = self.cost;
            nsbu = dims.nsbu;
            nsg = dims.nsg;
            if self.solver_options.N_horizon == 0
                return
            end

            nh_0 = dims.nh_0;
            nsh_0 = length(constraints.idxsh_0);
            if nsh_0 > nh_0
                error(['inconsistent dimension nsh_0 = ', num2str(nsh_0), '. Is greater than nh_0 = ', num2str(nh_0), '.']);
            end
            if any(constraints.idxsh_0 >= nh_0)
                error(['idxsh_0 = [', num2str(constraints.idxsh_0(:).'), '] contains value >= nh_0 = ', num2str(nh_0), '.']);
            end
            if any(constraints.idxsh_0 < 0)
                error(['idxsh_0 = [', num2str(constraints.idxsh_0(:).'), '] contains value < 0.']);
            end

            nphi_0 = dims.nphi_0;
            nsphi_0 = length(constraints.idxsphi_0);
            if nsphi_0 > nphi_0
                error(['inconsistent dimension nsphi_0 = ', num2str(nsphi_0), '. Is greater than nphi_0 = ', num2str(nphi_0), '.']);
            end
            if any(constraints.idxsphi_0 >= nphi_0)
                error(['idxsphi_0 = [', num2str(constraints.idxsphi_0(:).'), '] contains value >= nphi_0 = ', num2str(nphi_0), '.']);
            end
            if any(constraints.idxsphi_0 < 0)
                error(['idxsphi_0 = [', num2str(constraints.idxsphi_0(:).'), '] contains value < 0.']);
            end

            ns_0 = nsbu + nsg + nsh_0 + nsphi_0;
            wrong_field = '';
            if ns_0 == 0
                expected_shape = [0, 0];
            else
                expected_shape = [ns_0, 1];
            end

            if ~all(size(cost.Zl_0) == expected_shape)
                wrong_field = 'Zl_0';
                dim = size(cost.Zl_0);
            elseif ~all(size(cost.Zu_0) == expected_shape)
                wrong_field = 'Zu_0';
                dim = size(cost.Zu_0);
            elseif ~all(size(cost.zl_0) == expected_shape)
                wrong_field = 'zl_0';
                dim = size(cost.zl_0);
            elseif ~all(size(cost.zu_0) == expected_shape)
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
        end

        function make_consistent_slacks_terminal(self)
            constraints = self.constraints;
            dims = self.dims;
            cost = self.cost;

            nbx_e = dims.nbx_e;
            nsbx_e = length(constraints.idxsbx_e);
            if nsbx_e > nbx_e
                error(['inconsistent dimension nsbx_e = ', num2str(nsbx_e), '. Is greater than nbx_e = ', num2str(nbx_e), '.']);
            end
            if any(constraints.idxsbx_e >= nbx_e)
                error(['idxsbx_e = [', num2str(constraints.idxsbx_e(:).'), '] contains value >= nbx_e = ', num2str(nbx_e), '.']);
            end
            if any(constraints.idxsbx_e < 0)
                error(['idxsbx_e = [', num2str(constraints.idxsbx_e(:).'), '] contains value < 0.']);
            end

            ng_e = dims.ng_e;
            nsg_e = length(constraints.idxsg_e);
            if nsg_e > ng_e
                error(['inconsistent dimension nsg_e = ', num2str(nsg_e), '. Is greater than ng_e = ', num2str(ng_e), '.']);
            end
            if any(constraints.idxsg_e >= ng_e)
                error(['idxsg_e = [', num2str(constraints.idxsg_e(:).'), '] contains value >= ng_e = ', num2str(ng_e), '.']);
            end
            if any(constraints.idxsg_e < 0)
                error(['idxsg_e = [', num2str(constraints.idxsg_e(:).'), '] contains value < 0.']);
            end

            nh_e = dims.nh_e;
            nsh_e = length(constraints.idxsh_e);
            if nsh_e > nh_e
                error(['inconsistent dimension nsh_e = ', num2str(nsh_e), '. Is greater than nh_e = ', num2str(nh_e), '.']);
            end
            if any(constraints.idxsh_e >= nh_e)
                error(['idxsh_e = [', num2str(constraints.idxsh_e(:).'), '] contains value >= nh_e = ', num2str(nh_e), '.']);
            end
            if any(constraints.idxsh_e < 0)
                error(['idxsh_e = [', num2str(constraints.idxsh_e(:).'), '] contains value < 0.']);
            end


            nphi_e = dims.nphi_e;
            nsphi_e = length(constraints.idxsphi_e);
            if nsphi_e > nphi_e
                error(['inconsistent dimension nsphi_e = ', num2str(nsphi_e), '. Is greater than nphi_e = ', num2str(nphi_e), '.']);
            end
            if any(constraints.idxsphi_e >= nphi_e)
                error(['idxsphi_e = [', num2str(constraints.idxsphi_e(:).'), '] contains value >= nphi_e = ', num2str(nphi_e), '.']);
            end
            if any(constraints.idxsphi_e < 0)
                error(['idxsphi_e = [', num2str(constraints.idxsphi_e(:).'), '] contains value < 0.']);
            end


            ns_e = nsbx_e + nsg_e + nsh_e + nsphi_e;
            wrong_field = '';
            if ns_e == 0
                expected_shape = [0, 0];
            else
                expected_shape = [ns_e, 1];
            end

            if ~all(size(cost.Zl_e) == expected_shape)
                wrong_field = 'Zl_e';
                dim = size(cost.Zl_e);
            elseif ~all(size(cost.Zu_e) == expected_shape)
                wrong_field = 'Zu_e';
                dim = size(cost.Zu_e);
            elseif ~all(size(cost.zl_e) == expected_shape)
                wrong_field = 'zl_e';
                dim = size(cost.zl_e);
            elseif ~all(size(cost.zu_e) == expected_shape)
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

        end

        function make_consistent_discretization(self)
            dims = self.dims;
            opts = self.solver_options;

            if opts.N_horizon == 0
                opts.shooting_nodes = zeros(1, 1);
                opts.time_steps = ones(1, 1);
                return
            end

            if length(opts.tf) ~= 1 || opts.tf < 0
                error('time horizon tf should be a nonnegative number');
            end

            if ~isempty(opts.shooting_nodes)
                if opts.N_horizon + 1 ~= length(opts.shooting_nodes)
                    error('inconsistent dimension N regarding shooting nodes.');
                end
                for i=1:opts.N_horizon
                    opts.time_steps(i) = opts.shooting_nodes(i+1) - opts.shooting_nodes(i);
                end
                sum_time_steps = sum(opts.time_steps);
                if abs((sum_time_steps - opts.tf) / opts.tf) > 1e-14
                    warning('shooting nodes are not consistent with time horizon tf, rescaling automatically');
                    opts.time_steps = opts.time_steps * opts.tf / sum_time_steps;
                end
            elseif ~isempty(opts.time_steps)
                if opts.N_horizon ~= length(opts.time_steps)
                    error('inconsistent dimension N regarding time steps.');
                end
                sum_time_steps = sum(opts.time_steps);
                if abs((sum_time_steps - opts.tf) / opts.tf) > 1e-14
                    error(['time steps are not consistent with time horizon tf, ', ...
                        'got tf = ' num2str(opts.tf) '; sum(time_steps) = ' num2str(sum_time_steps) '.']);
                end
            else
                opts.time_steps = opts.tf/opts.N_horizon * ones(opts.N_horizon,1);
            end
            % add consistent shooting_nodes e.g. for plotting;
            if isempty(opts.shooting_nodes)
                opts.shooting_nodes = zeros(opts.N_horizon+1, 1);
                for i = 1:opts.N_horizon
                    opts.shooting_nodes(i+1) = sum(opts.time_steps(1:i));
                end
            end
            if any(opts.time_steps < 0)
                error(['ocp discretization: time_steps between shooting nodes must all be > 0', ...
                    ' got: ' num2str(opts.time_steps)])
            end
        end

        function make_consistent_simulation(self)
            opts = self.solver_options;
            if opts.N_horizon == 0
                return
            end

            % set integrator time automatically
            opts.Tsim = opts.time_steps(1);

            % integrator: num_stages
            if ~isempty(opts.sim_method_num_stages)
                if(strcmp(opts.integrator_type, "ERK"))
                    if (any(opts.sim_method_num_stages < 1) || any(opts.sim_method_num_stages > 4))
                        error(['ERK: num_stages = ', num2str(opts.sim_method_num_stages) ' not available. Only number of stages = {1,2,3,4} implemented!']);
                    end
                end
            end

            %% options sanity checks
            if length(opts.sim_method_num_steps) == 1
                opts.sim_method_num_steps = opts.sim_method_num_steps * ones(1, opts.N_horizon);
            elseif length(opts.sim_method_num_steps) ~= opts.N_horizon
                error('sim_method_num_steps must be a scalar or a vector of length N');
            end
            if length(opts.sim_method_num_stages) == 1
                opts.sim_method_num_stages = opts.sim_method_num_stages * ones(1, opts.N_horizon);
            elseif length(opts.sim_method_num_stages) ~= opts.N_horizon
                error('sim_method_num_stages must be a scalar or a vector of length N');
            end
            if length(opts.sim_method_jac_reuse) == 1
                opts.sim_method_jac_reuse = opts.sim_method_jac_reuse * ones(1, opts.N_horizon);
            elseif length(opts.sim_method_jac_reuse) ~= opts.N_horizon
                error('sim_method_jac_reuse must be a scalar or a vector of length N');
            end


        end

        function make_consistent(self, is_mocp_phase)
            if nargin < 2
                is_mocp_phase = false;
            end
            self.model.make_consistent(self.dims);

            model = self.model;
            dims = self.dims;
            cost = self.cost;
            constraints = self.constraints;
            opts = self.solver_options;

            self.detect_cost_and_constraints();

            if isempty(opts.N_horizon) && isempty(dims.N)
                error('N_horizon not provided.');
            elseif isempty(opts.N_horizon) && ~isempty(dims.N)
                opts.N_horizon = dims.N;
                disp(['field AcadosOcpDims.N has been migrated to AcadosOcpOptions.N_horizon.',...
                      ' setting AcadosOcpOptions.N_horizon = N.',...
                      ' For future comppatibility, please use AcadosOcpOptions.N_horizon directly.']);
            elseif ~isempty(opts.N_horizon) && ~isempty(dims.N) && opts.N_horizon ~= dims.N
                error(['Inconsistent dimension N, regarding N = ', num2str(dims.N),...
                       ', N_horizon = ', num2str(opts.N_horizon), '.']);
            else
                dims.N = opts.N_horizon;
            end

            % check if nx != nx_next
            if ~is_mocp_phase && dims.nx ~= dims.nx_next && opts.N_horizon > 1
                error(['nx_next = ', num2str(dims.nx_next), ' must be equal to nx = ', num2str(dims.nx), ' if more than one shooting interval is used.']);
            end

            % detect GNSF structure
            if strcmp(opts.integrator_type, 'GNSF') && opts.N_horizon > 0
                if dims.gnsf_nx1 + dims.gnsf_nx2 ~= dims.nx
                    % TODO: properly interface those.
                    gnsf_transcription_opts = struct();
                    detect_gnsf_structure(model, dims, gnsf_transcription_opts);
                else
                    warning('No GNSF model detected, assuming all required fields are set.')
                end
            end

            % sanity checks on options, which are done in setters in Python
            qp_solvers = {'PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES', 'FULL_CONDENSING_HPIPM', 'PARTIAL_CONDENSING_QPDUNES', 'PARTIAL_CONDENSING_OSQP', 'FULL_CONDENSING_DAQP'};
            if ~ismember(opts.qp_solver, qp_solvers)
                error(['Invalid qp_solver: ', opts.qp_solver, '. Available options are: ', strjoin(qp_solvers, ', ')]);
            end

            regularize_methods = {'NO_REGULARIZE', 'MIRROR', 'PROJECT', 'PROJECT_REDUC_HESS', 'CONVEXIFY', 'GERSHGORIN_LEVENBERG_MARQUARDT'};
            if ~ismember(opts.regularize_method, regularize_methods)
                error(['Invalid regularize_method: ', opts.regularize_method, '. Available options are: ', strjoin(regularize_methods, ', ')]);
            end
            hpipm_modes = {'BALANCE', 'SPEED_ABS', 'SPEED', 'ROBUST'};
            if ~ismember(opts.hpipm_mode, hpipm_modes)
                error(['Invalid hpipm_mode: ', opts.hpipm_mode, '. Available options are: ', strjoin(hpipm_modes, ', ')]);
            end
            INTEGRATOR_TYPES = {'ERK', 'IRK', 'GNSF', 'DISCRETE', 'LIFTED_IRK'};
            if ~ismember(opts.integrator_type, INTEGRATOR_TYPES)
                error(['Invalid integrator_type: ', opts.integrator_type, '. Available options are: ', strjoin(INTEGRATOR_TYPES, ', ')]);
            end

            COLLOCATION_TYPES = {'GAUSS_RADAU_IIA', 'GAUSS_LEGENDRE', 'EXPLICIT_RUNGE_KUTTA'};
            if ~ismember(opts.collocation_type, COLLOCATION_TYPES)
                error(['Invalid collocation_type: ', opts.collocation_type, '. Available options are: ', strjoin(COLLOCATION_TYPES, ', ')]);
            end

            COST_DISCRETIZATION_TYPES = {'EULER', 'INTEGRATOR'};
            if ~ismember(opts.cost_discretization, COST_DISCRETIZATION_TYPES)
                error(['Invalid cost_discretization: ', opts.cost_discretization, '. Available options are: ', strjoin(COST_DISCRETIZATION_TYPES, ', ')]);
            end

            search_direction_modes = {'NOMINAL_QP', 'BYRD_OMOJOKUN', 'FEASIBILITY_QP'};
            if ~ismember(opts.search_direction_mode, search_direction_modes)
                error(['Invalid search_direction_mode: ', opts.search_direction_mode, '. Available options are: ', strjoin(search_direction_modes, ', ')]);
            end

            qpscaling_scale_constraints_types = {'INF_NORM', 'NO_CONSTRAINT_SCALING'};
            if ~ismember(opts.qpscaling_scale_constraints, qpscaling_scale_constraints_types)
                error(['Invalid qpscaling_scale_constraints: ', opts.qpscaling_scale_constraints, '. Available options are: ', strjoin(qpscaling_scale_constraints_types, ', ')]);
            end

            qpscaling_scale_objective_types = {'OBJECTIVE_GERSHGORIN', 'NO_OBJECTIVE_SCALING'};
            if ~ismember(opts.qpscaling_scale_objective, qpscaling_scale_objective_types)
                error(['Invalid qpscaling_scale_objective: ', opts.qpscaling_scale_objective, '. Available options are: ', strjoin(qpscaling_scale_objective_types, ', ')]);
            end
            % OCP name
            self.name = model.name;

            % parameters
            if isempty(self.parameter_values)
                if dims.np > 0
                    warning(['self.parameter_values are not set.', ...
                            10 'Using zeros(np,1) by default.' 10 'You can update them later using set().']);
                end
                self.parameter_values = zeros(self.dims.np,1);
            elseif length(self.parameter_values) ~= self.dims.np
                error(['parameter_values has the wrong shape. Expected: ' num2str(self.dims.np)])
            end


            % parameters
            if isempty(self.p_global_values)
                if dims.np_global > 0
                    warning(['self.p_global_values are not set.', ...
                            10 'Using zeros(np_global,1) by default.' 10 'You can update them later using set().']);
                end
                self.p_global_values = zeros(self.dims.np_global,1);
            elseif length(self.p_global_values) ~= self.dims.np_global
                error(['p_global_values has the wrong shape. Expected: ' num2str(self.dims.np_global)])
            end

            %% cost
            self.make_consistent_cost_initial();
            self.make_consistent_cost_path();
            self.make_consistent_cost_terminal();

            % cost integration
            if strcmp(opts.cost_discretization, "INTEGRATOR") && opts.N_horizon > 0
                if ~(strcmp(cost.cost_type, "NONLINEAR_LS") || strcmp(cost.cost_type, "CONVEX_OVER_NONLINEAR"))
                    error('INTEGRATOR cost discretization requires CONVEX_OVER_NONLINEAR or NONLINEAR_LS cost type for path cost.')
                end
                if ~(strcmp(cost.cost_type_0, "NONLINEAR_LS") || strcmp(cost.cost_type_0, "CONVEX_OVER_NONLINEAR"))
                    error('INTEGRATOR cost discretization requires CONVEX_OVER_NONLINEAR or NONLINEAR_LS cost type for initial cost.')
                end
                if strcmp(opts.nlp_solver_type, 'SQP_WITH_FEASIBLE_QP')
                    error('cost_discretization == INTEGRATOR is not compatible with SQP_WITH_FEASIBLE_QP yet.')
                end
            end


            %% constraints
            % qpdunes
            if ~isempty(strfind(opts.qp_solver, 'QPDUNES'))
                constraints.idxbxe_0 = [];
                dims.nbxe_0 = 0;
            end
            self.make_consistent_constraints_initial();
            self.make_consistent_constraints_path();
            self.make_consistent_constraints_terminal();

            %% slack dimensions
            self.make_consistent_slacks_path();
            self.make_consistent_slacks_initial();
            self.make_consistent_slacks_terminal();

            % check for ACADOS_INFTY
            if ~ismember(opts.qp_solver, {'PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_DAQP'})
                ACADOS_INFTY = get_acados_infty();
                % loop over all bound vectors
                if opts.N_horizon > 0
                    fields = {'lbx_e', 'ubx_e', 'lg_e', 'ug_e', 'lh_e', 'uh_e', 'lphi_e', 'uphi_e'};
                else
                    fields = {'lbx_0', 'ubx_0', 'lbx', 'ubx', 'lbx_e', 'ubx_e', 'lg', 'ug', 'lg_e', 'ug_e', 'lh', 'uh', 'lh_e', 'uh_e', 'lbu', 'ubu', 'lphi', 'uphi', 'lphi_e', 'uphi_e'};
                end
                for i = 1:length(fields)
                    field = fields{i};
                    bound = constraints.(field);
                    if any(bound >= ACADOS_INFTY) || any(bound <= -ACADOS_INFTY)
                        error(['Field ', field, ' contains values outside the interval (-ACADOS_INFTY, ACADOS_INFTY) with ACADOS_INFTY = ', num2str(ACADOS_INFTY, '%.2e'), '. One-sided constraints are not supported by the chosen QP solver ', opts.qp_solver, '.']);
                    end
                end
            end

            self.make_consistent_discretization();

            % cost_scaling
            if isempty(opts.cost_scaling)
                opts.cost_scaling = [opts.time_steps(:); 1.0];
            elseif length(opts.cost_scaling) ~= opts.N_horizon+1
                error(['cost_scaling must have length N+1 = ', num2str(N+1)]);
            end

            self.make_consistent_simulation();

            if strcmp(opts.qp_solver, "PARTIAL_CONDENSING_HPMPC") || ...
                strcmp(opts.qp_solver, "PARTIAL_CONDENSING_QPDUNES") || ...
                strcmp(opts.qp_solver, "PARTIAL_CONDENSING_OOQP")
                if self.dims.ns > 0 || self.dims.ns_e > 0
                    error(['selected QP solver ', opts.qp_solver, ' does not support soft constraints (yet).'])
                end
            end

            if ~(strcmp(opts.qp_solver, "FULL_CONDENSING_HPIPM") || ...
                strcmp(opts.qp_solver, "PARTIAL_CONDENSING_HPIPM") || ...
                strcmp(opts.qp_solver, "FULL_CONDENSING_DAQP"))
                disp(['NOTE: The selected QP solver ', opts.qp_solver, ' does not support one-sided constraints yet.']);
            end

            % fixed hessian
            if opts.fixed_hess
                if opts.hessian_approx == 'EXACT'
                    error('fixed_hess and hessian_approx = EXACT are incompatible')
                end
                if ~strcmp(cost.cost_type_0, "LINEAR_LS") && opts.N_horizon > 0
                    error('fixed_hess requires LINEAR_LS cost_type_0')
                end
                if ~strcmp(cost.cost_type, "LINEAR_LS") && opts.N_horizon > 0
                    error('fixed_hess requires LINEAR_LS cost_type')
                end
                if ~strcmp(cost.cost_type_e, "LINEAR_LS")
                    error('fixed_hess requires LINEAR_LS cost_type_e')
                end
            end

            % TODO: add checks for solution sensitivities when brining them to MATLAB

            % check if qp_solver_cond_N is set
            if isempty(opts.qp_solver_cond_N)
                opts.qp_solver_cond_N = opts.N_horizon;
            end
            if opts.qp_solver_cond_N > opts.N_horizon
                error('qp_solver_cond_N > N_horizon is not supported.');
            end

            if ~isempty(opts.qp_solver_cond_block_size)
                if sum(opts.qp_solver_cond_block_size) ~= opts.N_horizon
                    error(['sum(qp_solver_cond_block_size) =', num2str(sum(opts.qp_solver_cond_block_size)), ' != N = {opts.N_horizon}.']);
                end
                if length(opts.qp_solver_cond_block_size) ~= opts.qp_solver_cond_N+1
                    error('qp_solver_cond_block_size should have length qp_solver_cond_N+1.');
                end
            end

            if strcmp(opts.nlp_solver_type, "DDP")
                if opts.N_horizon == 0
                    error('DDP solver only supported for N_horizon > 0.');
                end
                if ~strcmp(opts.qp_solver, "PARTIAL_CONDENSING_HPIPM") || (opts.qp_solver_cond_N ~= opts.N_horizon)
                    error('DDP solver only supported for PARTIAL_CONDENSING_HPIPM with qp_solver_cond_N == N_horizon.');
                end
                if any([dims.nbu, dims.nbx, dims.ng, dims.nh, dims.nphi])
                    error('DDP only supports initial state constraints, got path constraints.')
                end
                if any([dims.ng_e, dims.nphi_e, dims.nh_e])
                    error('DDP only supports initial state constraints, got terminal constraints.')
                end
            end

            if ~ismember(opts.qp_solver_t0_init, [0, 1, 2])
                error('qp_solver_t0_init must be one of [0, 1, 2].');
            end

            if opts.tau_min > 0 && isempty(strfind(opts.qp_solver, 'HPIPM'))
                error('tau_min > 0 is only compatible with HPIPM.');
            end

            if opts.N_horizon == 0
                cost_types_to_check = [strcmp(cost.cost_type_e, {'LINEAR_LS', 'NONLINEAR_LS'})];
            else
                cost_types_to_check = [strcmp(cost.cost_type, {'LINEAR_LS', 'NONLINEAR_LS'}) ...
                                            strcmp(cost.cost_type_0, {'LINEAR_LS', 'NONLINEAR_LS'}) ...
                                            strcmp(cost.cost_type_e, {'LINEAR_LS', 'NONLINEAR_LS'})];
            end
            if (opts.as_rti_level == 1 || opts.as_rti_level == 2) && any(cost_types_to_check)
                error('as_rti_level in [1, 2] not supported for LINEAR_LS and NONLINEAR_LS cost type.');
            end

            if ~strcmp(opts.qpscaling_scale_constraints, "NO_CONSTRAINT_SCALING") || ~strcmp(opts.qpscaling_scale_objective, "NO_OBJECTIVE_SCALING")
                if strcmp(opts.nlp_solver_type, "SQP_RTI")
                    error('qpscaling_scale_constraints and qpscaling_scale_objective not supported for SQP_RTI solver.');
                end
            end

            % Set default parameters for globalization
            ddp_with_merit_or_funnel = strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH') || (strcmp(opts.globalization, 'MERIT_BACKTRACKING') && strcmp(opts.nlp_solver_type, 'DDP'));

            if isempty(opts.globalization_alpha_min)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_alpha_min = 1e-17;
                else
                    opts.globalization_alpha_min = 0.05;
                end
            end

            if isempty(opts.globalization_alpha_reduction)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_alpha_reduction = 0.5;
                else
                    opts.globalization_alpha_reduction = 0.7;
                end
            end

            if isempty(opts.globalization_eps_sufficient_descent)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_eps_sufficient_descent = 1e-6;
                else
                    opts.globalization_eps_sufficient_descent = 1e-4;
                end
            end

            if isempty(opts.eval_residual_at_max_iter)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.eval_residual_at_max_iter = true;
                else
                    opts.eval_residual_at_max_iter = false;
                end
            end

            if isempty(opts.globalization_full_step_dual)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_full_step_dual = 1;
                else
                    opts.globalization_full_step_dual = 0;
                end
            end

            % sanity check for Funnel globalization and SQP
            if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH') && ~strcmp(opts.nlp_solver_type, 'SQP')
                error('FUNNEL_L1PEN_LINESEARCH only supports SQP.');
            end

            % termination
            if isempty(opts.nlp_solver_tol_min_step_norm)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.nlp_solver_tol_min_step_norm = 1e-12;
                else
                    opts.nlp_solver_tol_min_step_norm = 0.0;
                end
            end

            %% Deprecated / migrated options
            if ~isempty(opts.nlp_solver_step_length)
                warning('nlp_solver_step_length is deprecated, use globalization_fixed_step_length instead.');
                if opts.globalization_fixed_step_length ~= 1.0
                    error('nlp_solver_step_length and globalization_fixed_step_length are both set, please use only globalization_fixed_step_length.');
                end
                opts.globalization_fixed_step_length = opts.nlp_solver_step_length;
            end

            if opts.globalization_fixed_step_length < 0.0 || opts.globalization_fixed_step_length > 1.0
                error('globalization_fixed_step_length must be in [0, 1].');
            end

            % Set default parameters for globalization
            if isempty(opts.globalization_alpha_min)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_alpha_min = 1e-17;
                else
                    opts.globalization_alpha_min = 0.05;
                end
            end

            if isempty(opts.globalization_alpha_reduction)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_alpha_reduction = 0.5;
                else
                    opts.globalization_alpha_reduction = 0.7;
                end
            end

            if isempty(opts.globalization_eps_sufficient_descent)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_eps_sufficient_descent = 1e-6;
                else
                    opts.globalization_eps_sufficient_descent = 1e-4;
                end
            end

            if isempty(opts.globalization_full_step_dual)
                % if strcmp(opts.globalization, 'FUNNEL_L1PEN_LINESEARCH')
                if ddp_with_merit_or_funnel
                    opts.globalization_full_step_dual = 1;
                else
                    opts.globalization_full_step_dual = 0;
                end
            end

            if isa(self.zoro_description, 'ZoroDescription')
                if opts.N_horizon == 0
                    error('ZORO only supported for N_horizon > 0.');
                end
                self.zoro_description.process();
            end

            % Anderson acceleration
            if opts.with_anderson_acceleration
                if strcmp(opts.nlp_solver_type, "DDP")
                    error('Anderson acceleration not supported for DDP solver.');
                end
                if ~strcmp(opts.globalization, "FIXED_STEP")
                    error('Anderson acceleration only supported for FIXED_STEP globalization for now.');
                end
            end

            % check terminal stage
            fields = {'cost_expr_ext_cost_e', 'cost_expr_ext_cost_custom_hess_e', ...
                      'cost_y_expr_e', 'cost_psi_expr_e', 'cost_conl_custom_outer_hess_e', ...
                      'con_h_expr_e', 'con_phi_expr_e', 'con_r_expr_e'};
            for i = 1:length(fields)
                field = fields{i};
                val = model.(field);
                if ~isempty(val) && (depends_on(val, model.u) || depends_on(val, model.z))
                    error([field ' can not depend on u or z.'])
                end
            end
        end

        function [] = detect_cost_and_constraints(self)
            % detect cost type
            N = self.solver_options.N_horizon;
            if N == 0
                if strcmp(self.cost.cost_type_e, 'AUTO')
                    detect_cost_type(self.model, self.cost, self.dims, 'terminal');
                end
                if strcmp(self.constraints.constr_type_e, 'AUTO')
                    detect_constraint_structure(self.model, self.constraints, 'terminal');
                end
                return
            end

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
                self.cost.Vx_0 = self.cost.Vx;
                self.cost.Vu_0 = self.cost.Vu;
                self.cost.Vz_0 = self.cost.Vz;
                self.model.cost_y_expr_0 = self.model.cost_y_expr;
                self.cost.cost_ext_fun_type_0 = self.cost.cost_ext_fun_type;
                self.model.cost_expr_ext_cost_0 = self.model.cost_expr_ext_cost;
                self.model.cost_expr_ext_cost_custom_hess_0 = self.model.cost_expr_ext_cost_custom_hess;
                self.cost.cost_source_ext_cost_0 = self.cost.cost_source_ext_cost;
                self.cost.cost_function_ext_cost_0 = self.cost.cost_function_ext_cost;
                self.cost.W_0 = self.cost.W;
                self.cost.yref_0 = self.cost.yref;
                self.model.cost_psi_expr_0 = self.model.cost_psi_expr;
                self.model.cost_r_in_psi_expr_0 = self.model.cost_r_in_psi_expr;
            end

            % detect constraint structure
            constraint_types = {self.constraints.constr_type_0, self.constraints.constr_type, self.constraints.constr_type_e};
            for n=1:3
                if strcmp(constraint_types{n}, 'AUTO')
                    detect_constraint_structure(self.model, self.constraints, stage_types{n});
                end
            end
        end

        function context = generate_external_functions(ocp, context)

            %% generate C code for CasADi functions / copy external functions
            solver_opts = ocp.solver_options;

            if nargin < 2
                % options for code generation
                code_gen_opts = struct();
                code_gen_opts.generate_hess = strcmp(solver_opts.hessian_approx, 'EXACT');
                code_gen_opts.with_solution_sens_wrt_params = solver_opts.with_solution_sens_wrt_params;
                code_gen_opts.with_value_sens_wrt_params = solver_opts.with_value_sens_wrt_params;
                code_gen_opts.code_export_directory = ocp.code_export_directory;

                code_gen_opts.ext_fun_expand_dyn = solver_opts.ext_fun_expand_dyn;
                code_gen_opts.ext_fun_expand_cost = solver_opts.ext_fun_expand_cost;
                code_gen_opts.ext_fun_expand_constr = solver_opts.ext_fun_expand_constr;
                code_gen_opts.ext_fun_expand_precompute = solver_opts.ext_fun_expand_precompute;

                context = GenerateContext(ocp.model.p_global, ocp.name, code_gen_opts);
            else
                code_gen_opts = context.opts;
            end
            context = setup_code_generation_context(ocp, context, false, false);
            context.finalize();
            ocp.external_function_files_model = context.get_external_function_file_list(false);
            ocp.external_function_files_ocp = context.get_external_function_file_list(true);
            ocp.dims.n_global_data = context.get_n_global_data();
        end

        function context = setup_code_generation_context(ocp, context, ignore_initial, ignore_terminal)
            code_gen_opts = context.opts;
            solver_opts = ocp.solver_options;
            constraints = ocp.constraints;
            cost = ocp.cost;
            dims = ocp.dims;

            setup_code_generation_context_dynamics(ocp, context);

            if solver_opts.N_horizon == 0
                stage_type_indices = [3];
            else
                if ignore_initial && ignore_terminal
                    stage_type_indices = [2];
                elseif ignore_terminal
                    stage_type_indices = [1, 2];
                elseif ignore_initial
                    stage_type_indices = [2, 3];
                else
                    stage_type_indices = [1, 2, 3];
                end
            end

            stage_types = {'initial', 'path', 'terminal'};

            % cost
            cost_types = {cost.cost_type_0, cost.cost_type, cost.cost_type_e};
            cost_ext_fun_types = {cost.cost_ext_fun_type_0, cost.cost_ext_fun_type, cost.cost_ext_fun_type_e};
            cost_dir = fullfile(pwd, ocp.code_export_directory, [ocp.name '_cost']);

            for n = 1:length(stage_type_indices)

                i = stage_type_indices(n);
                if strcmp(cost_ext_fun_types{i}, 'generic')
                    if strcmp(cost_types{i}, 'EXTERNAL')
                        setup_generic_cost(context, cost, cost_dir, stage_types{i})
                    else
                        error('Unknown cost_type for cost_ext_fun_types generic: got %s', cost_types{i});
                    end

                else
                    check_casadi_version();
                    switch cost_types{i}
                        case 'LINEAR_LS'
                            continue
                        case 'NONLINEAR_LS'
                            generate_c_code_nonlinear_least_squares(context, ocp.model, cost_dir, stage_types{i});

                        case 'CONVEX_OVER_NONLINEAR'
                            error("Convex-over-nonlinear cost is not implemented yet.")

                        case 'EXTERNAL'
                            generate_c_code_ext_cost(context, ocp.model, cost_dir, stage_types{i});
                        otherwise
                            error('Unknown value for cost_ext_fun_types %s', cost_ext_fun_types{i});
                    end
                end
            end

            % constraints
            constraints_types = {constraints.constr_type_0, constraints.constr_type, constraints.constr_type_e};
            constraints_dims = {dims.nh_0, dims.nh, dims.nh_e};
            constraints_dir = fullfile(pwd, ocp.code_export_directory, [ocp.name '_constraints']);

            for n = 1:length(stage_type_indices)
                i = stage_type_indices(n);
                if strcmp(constraints_types{i}, 'BGH') && constraints_dims{i} > 0
                    generate_c_code_nonlinear_constr(context, ocp.model, constraints_dir, stage_types{i});
                end
            end
        end

        function setup_code_generation_context_dynamics(ocp, context)
            code_gen_opts = context.opts;
            solver_opts = ocp.solver_options;
            if solver_opts.N_horizon == 0
                return
            end

            model_dir = fullfile(pwd, code_gen_opts.code_export_directory, [ocp.name '_model']);

            if strcmp(ocp.model.dyn_ext_fun_type, 'generic')
                check_dir_and_create(model_dir);
                copyfile(fullfile(pwd, ocp.model.dyn_generic_source), model_dir);
                context.add_external_function_file(ocp.model.dyn_generic_source, model_dir);
            elseif strcmp(ocp.model.dyn_ext_fun_type, 'casadi')
                check_casadi_version();
                switch solver_opts.integrator_type
                    case 'ERK'
                        generate_c_code_explicit_ode(context, ocp.model, model_dir);
                    case 'IRK'
                        generate_c_code_implicit_ode(context, ocp.model, model_dir);
                    case 'LIFTED_IRK'
                        if ~(isempty(ocp.model.t) || length(ocp.model.t) == 0)
                            error('NOT LIFTED_IRK with time-varying dynamics not implemented yet.')
                        end
                        generate_c_code_implicit_ode(context, ocp.model, model_dir);
                    case 'GNSF'
                        generate_c_code_gnsf(context, ocp.model, model_dir);
                    case 'DISCRETE'
                        generate_c_code_discrete_dynamics(context, ocp.model, model_dir);
                    otherwise
                        error('Unknown integrator type.')
                end
            else
                error('Unknown dyn_ext_fun_type.')
            end
        end

        function render_templates(self)

            %% render templates
            json_fullfile = fullfile(pwd, self.json_file);
            main_dir = pwd;
            chdir(self.code_export_directory);

            template_list = self.get_template_list();
            for i = 1:length(template_list)
                in_file = template_list{i}{1};
                out_file = template_list{i}{2};
                if length(template_list{i}) == 3
                    out_dir = template_list{i}{3};
                    if ~(exist(out_dir, 'dir'))
                        mkdir(out_dir);
                    end
                    out_file = fullfile(out_dir, out_file);
                end
                render_file( in_file, out_file, json_fullfile )
            end

            fprintf('Successfully rendered acados templates!\n');


            % Custom templates
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            custom_template_glob = fullfile(acados_folder, 'interfaces', 'acados_template', 'acados_template', 'custom_update_templates', '*');

            template_list = self.solver_options.custom_templates;
            for i = 1:length(self.solver_options.custom_templates)
                in_file = template_list{i}{1};
                out_file = template_list{i}{2};
                if length(template_list{i}) == 3
                    out_dir = template_list{i}{3};
                    if ~(exist(out_dir, 'dir'))
                        mkdir(out_dir);
                    end
                    out_file = fullfile(out_dir, out_file);
                end
                render_file( in_file, out_file, json_fullfile, custom_template_glob )
            end

            cd(main_dir)
        end


        function template_list = get_template_list(self)
            % returns a cell of cells in the form:
            % (input_filename, output_filname)
            % or
            % (input_filename, output_filname, output_directory)
            template_list = {};
            template_list{end+1} = {'main.in.c', ['main_', self.name, '.c']};
            template_list{end+1} = {'acados_solver.in.h', ['acados_solver_', self.name, '.h']};
            template_list{end+1} = {'acados_solver.in.c', ['acados_solver_', self.name, '.c']};
            template_list{end+1} = {'CMakeLists.in.txt', ['CMakeLists.txt']};
            template_list{end+1} = {'Makefile.in', ['Makefile']};

            % integrator
            if ~strcmp(self.solver_options.integrator_type, 'DISCRETE')
                template_list{end+1} = {'acados_sim_solver.in.c', ['acados_sim_solver_', self.name, '.c']};
                template_list{end+1} = {'acados_sim_solver.in.h', ['acados_sim_solver_', self.name, '.h']};
                template_list{end+1} = {'main_sim.in.c', ['main_sim_', self.name, '.c']};
            end
            % MEX files
            matlab_template_path = 'matlab_templates';
            template_list{end+1} = {fullfile(matlab_template_path, 'mex_solver.in.m'), [self.name, '_mex_solver.m']};
            template_list{end+1} = {fullfile(matlab_template_path, 'make_mex.in.m'), ['make_mex_', self.name, '.m']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_create.in.c'), ['acados_mex_create_', self.name, '.c']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_free.in.c'), ['acados_mex_free_', self.name, '.c']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_solve.in.c'), ['acados_mex_solve_', self.name, '.c']};
            template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_set.in.c'), ['acados_mex_set_', self.name, '.c']};

            if ~isempty(self.solver_options.custom_update_filename)
                template_list{end+1} = {fullfile(matlab_template_path, 'acados_mex_custom_update.in.c'), ['acados_mex_custom_update_', self.name, '.c']};
            end

            % append headers
            template_list = [template_list, self.get_external_function_header_templates()];

            if self.dims.n_global_data > 0
                template_list{end+1} = {'p_global_precompute_fun.in.h',  [self.model.name, '_p_global_precompute_fun.h']};
            end

            % Simulink
            if ~isempty(self.simulink_opts)
                template_list{end+1} = {fullfile(matlab_template_path, 'acados_solver_sfun.in.c'), ['acados_solver_sfunction_', self.name, '.c']};
                template_list{end+1} = {fullfile(matlab_template_path, 'make_sfun.in.m'), ['make_sfun.m']};
                if ~strcmp(self.solver_options.integrator_type, 'DISCRETE')
                    template_list{end+1} = {fullfile(matlab_template_path, 'acados_sim_solver_sfun.in.c'), ['acados_sim_solver_sfunction_', self.name, '.c']};
                    template_list{end+1} = {fullfile(matlab_template_path, 'make_sfun_sim.in.m'), ['make_sfun_sim.m']};
                end
                if self.simulink_opts.inputs.rti_phase && ~strcmp(self.solver_options.nlp_solver_type, 'SQP_RTI')
                    error('rti_phase is only supported for SQP_RTI');
                end
                if self.simulink_opts.outputs.KKT_residuals && strcmp(self.solver_options.nlp_solver_type, 'SQP_RTI')
                    warning('KKT_residuals now computes the residuals of the output iterate in SQP_RTI, this leads to increased computation time, turn off this port if it is not needed. See https://github.com/acados/acados/pull/1346.');
                end
            else
                disp("Not rendering Simulink-related templates, as simulink_opts are not specified.")
            end
        end

        function template_list = get_external_function_header_templates(self)
            template_list = {};
            template_list{end+1} = {'model.in.h', [self.model.name, '_model.h'], [self.model.name, '_model']};
            template_list{end+1} = {'cost.in.h', [self.model.name, '_cost.h'], [self.model.name, '_cost']};
            template_list{end+1} = {'constraints.in.h', [self.model.name, '_constraints.h'], [self.model.name, '_constraints']};
        end

        function dump_to_json(self, json_file)
            if nargin < 2
                json_file = self.json_file;
            end

            out_struct = orderfields(self.struct());

            % add compilation information to json
            acados_folder = getenv('ACADOS_INSTALL_DIR');
            libs = loadjson(fileread(fullfile(acados_folder, 'lib', 'link_libs.json')));
            out_struct.acados_link_libs = orderfields(libs);
            if ismac
                out_struct.os = 'mac';
            elseif isunix
                out_struct.os = 'unix';
            else
                out_struct.os = 'pc';
            end

            % prepare struct for json dump
            out_struct.parameter_values = reshape(num2cell(self.parameter_values), [1, self.dims.np]);
            out_struct.p_global_values = reshape(num2cell(self.p_global_values), [1, self.dims.np_global]);
            out_struct.model = orderfields(self.model.convert_to_struct_for_json_dump());
            out_struct.dims = orderfields(out_struct.dims.struct());
            out_struct.cost = orderfields(out_struct.cost.convert_to_struct_for_json_dump());
            out_struct.constraints = orderfields(out_struct.constraints.convert_to_struct_for_json_dump());
            out_struct.solver_options = orderfields(out_struct.solver_options.convert_to_struct_for_json_dump(self.solver_options.N_horizon));

            if ~isempty(self.zoro_description)
                out_struct.zoro_description = orderfields(self.zoro_description.convert_to_struct_for_json_dump());
            end

            % actual json dump
            json_string = savejson('', out_struct, 'ForceRootName', 0);
            fid = fopen(json_file, 'w');
            if fid == -1, error('Cannot create JSON file'); end
            fwrite(fid, json_string, 'char');
            fclose(fid);
        end
    end % methods
end


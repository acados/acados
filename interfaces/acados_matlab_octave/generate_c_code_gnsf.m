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


function generate_c_code_gnsf(context, model, model_dir)

    import casadi.*

    A = model.dyn_gnsf_A;
    B = model.dyn_gnsf_B;
    C = model.dyn_gnsf_C;
    E = model.dyn_gnsf_E;
    c = model.dyn_gnsf_c;

    L_x    = model.dyn_gnsf_L_x;
    L_z    = model.dyn_gnsf_L_z;
    L_xdot = model.dyn_gnsf_L_xdot;
    L_u    = model.dyn_gnsf_L_u;

    A_LO = model.dyn_gnsf_A_LO;
    E_LO = model.dyn_gnsf_E_LO;
    B_LO = model.dyn_gnsf_B_LO;
    c_LO = model.dyn_gnsf_c_LO;

    % state permutation vector: x_gnsf = dvecpe(x, ipiv)
    ipiv_x = model.dyn_gnsf_ipiv_x;
    idx_perm_x = model.dyn_gnsf_idx_perm_x;
    ipiv_z = model.dyn_gnsf_ipiv_z;
    idx_perm_z = model.dyn_gnsf_idx_perm_z;
    ipiv_f = model.dyn_gnsf_ipiv_f;
    idx_perm_f = model.dyn_gnsf_idx_perm_f;

    % expressions
    phi = model.dyn_gnsf_expr_phi;
    f_lo = model.dyn_gnsf_expr_f_lo;

    % binaries
    nontrivial_f_LO = model.gnsf_nontrivial_f_LO;
    purely_linear = model.gnsf_purely_linear;

    % symbolics
    y = model.sym_gnsf_y;
    uhat = model.sym_gnsf_uhat;
    x = model.x;
    xdot = model.xdot;
    u = model.u;
    z = model.z;
    p = model.p;

    model_name = model.name;

    nx1 = size(L_x, 2);
    nz1 = size(L_z, 2);

    % CasADi variables and expressions
    x1 = x(idx_perm_x(1:nx1));
    x1dot = xdot(idx_perm_x(1:nx1));
    z1 = z(idx_perm_z(1:nz1));

    %% generate functions
    if ~purely_linear
        jac_phi_y = jacobian(phi,y);
        jac_phi_uhat = jacobian(phi,uhat);

        context.add_function_definition([model_name,'_gnsf_phi_fun'], {y, uhat, p}, {phi}, model_dir, 'dyn');
        context.add_function_definition([model_name,'_gnsf_phi_fun_jac_y'], {y, uhat, p}, {phi, jac_phi_y}, model_dir, 'dyn');
        context.add_function_definition([model_name,'_gnsf_phi_jac_y_uhat'], {y, uhat, p}, {jac_phi_y, jac_phi_uhat}, model_dir, 'dyn');


        if nontrivial_f_LO
            context.add_function_definition([model_name,'_gnsf_f_lo_fun_jac_x1k1uz'], {x1, x1dot, z1, u, p}, ...
                {f_lo, [jacobian(f_lo,x1), jacobian(f_lo,x1dot), jacobian(f_lo,u), jacobian(f_lo,z1)]}, model_dir, 'dyn');
        end
    end

    % get_matrices function
    dummy = x(1);
    context.add_function_definition([model_name,'_gnsf_get_matrices_fun'], {dummy},...
        {A, B, C, E, L_x, L_xdot, L_z, L_u, A_LO, c, E_LO, B_LO,...
        nontrivial_f_LO, purely_linear, ipiv_x, ipiv_z, c_LO}, model_dir, 'dyn');
end

%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%


function dump_gnsf_functions(model)

    acados_folder = getenv('ACADOS_INSTALL_DIR');
    addpath(fullfile(acados_folder, 'external', 'jsonlab'))

    %% import casadi
    import casadi.*

    casadi_version = CasadiMeta.version();
    if ~(strcmp(casadi_version(1:3),'3.5') || strcmp(casadi_version(1:3),'3.6'))
        warning('CasADi serialization requires CasADi >= 3.5, you are using: %s.', casadi_version);
    end

    %% import models
    % model matrices
    A  = model.gnsf_model.A;
    B  = model.gnsf_model.B;
    C  = model.gnsf_model.C;
    E  = model.gnsf_model.E;
    c  = model.gnsf_model.c;

    L_x    = model.gnsf_model.L_x;
    L_z    = model.gnsf_model.L_z;
    L_xdot = model.gnsf_model.L_xdot;
    L_u    = model.gnsf_model.L_u;

    A_LO = model.gnsf_model.A_LO;
    E_LO = model.gnsf_model.E_LO;
    B_LO = model.gnsf_model.B_LO;
    c_LO = model.gnsf_model.c_LO;

    % state permutation vector: x_gnsf = dvecpe(x, ipiv)
    ipiv_x = model.gnsf_model.ipiv_x;
    idx_perm_x = model.gnsf_model.idx_perm_x;
    ipiv_z = model.gnsf_model.ipiv_z;
    idx_perm_z = model.gnsf_model.idx_perm_z;
    ipiv_f = model.gnsf_model.ipiv_f;
    idx_perm_f = model.gnsf_model.idx_perm_f;

    % CasADi variables and expressions
    % x
    x = model.x;
    x1 = x(model.gnsf_model.idx_perm_x(1:model.gnsf_model.dims.nx1));
    % check type
    if isa(x(1), 'casadi.SX')
        isSX = true;
    else
        isSX = false;
    end
    % xdot
    xdot = model.xdot;
    x1dot = xdot(model.gnsf_model.idx_perm_x(1:model.gnsf_model.dims.nx1));
    u = model.u;
    if length(model.z) > 0
        z = model.z;
        z1 = model.z(model.gnsf_model.idx_perm_z(1:model.gnsf_model.dims.nz1));
    else
        if isSX
            z = SX.sym('z',0, 0);
            z1 = SX.sym('z1',0, 0);
        else
            z = MX.sym('z',0, 0);
            z1 = MX.sym('z1',0, 0);
        end
    end

    p = model.p;
    y = model.sym_gnsf_y;
    uhat = model.sym_gnsf_uhat;

    % expressions
    phi = model.gnsf_model.phi;
    f_lo = model.gnsf_model.f_lo;

    nontrivial_f_LO = model.model_gnsf.nontrivial_f_LO;
    purely_linear = model.model_gnsf.purely_linear;

    %% generate functions
    if ~purely_linear
        jac_phi_y = jacobian(phi,y);
        jac_phi_uhat = jacobian(phi,uhat);

        phi_fun = Function([model.name,'_gnsf_phi_fun'], {y, uhat, p}, {phi});
        phi_fun_jac_y = Function([model.name,'_gnsf_phi_fun_jac_y'], {y, uhat, p}, {phi, jac_phi_y});
        phi_jac_y_uhat = Function([model.name,'_gnsf_phi_jac_y_uhat'], {y, uhat, p}, {jac_phi_y, jac_phi_uhat});

        if nontrivial_f_LO
            f_lo_fun_jac_x1k1uz = Function([model.name,'_gnsf_f_lo_fun_jac_x1k1uz'], {x1, x1dot, z1, u, p}, ...
                {f_lo, [jacobian(f_lo,x1), jacobian(f_lo,x1dot), jacobian(f_lo,u), jacobian(f_lo,z1)]});
        end
    end

    % get_matrices function
    dummy = x(1);
    get_matrices_fun = Function([model.name,'_gnsf_get_matrices_fun'], {dummy},...
        {A, B, C, E, L_x, L_xdot, L_z, L_u, A_LO, c, E_LO, B_LO,...
        nontrivial_f_LO, purely_linear, ipiv_x, ipiv_z, c_LO});

    % dump functions to json
    out = struct();
    out.phi_fun = phi_fun.serialize();
    out.phi_fun_jac_y = phi_fun_jac_y.serialize();
    out.phi_jac_y_uhat = phi_jac_y_uhat.serialize();
    if exist('f_lo_fun_jac_x1k1uz', 'var')
        out.f_lo_fun_jac_x1k1uz = f_lo_fun_jac_x1k1uz.serialize();
    end
    out.get_matrices_fun = get_matrices_fun.serialize();
    out.casadi_version = casadi_version;

    json_filename = [model.name '_gnsf_functions.json'];
    json_string = savejson('', out, 'ForceRootName', 0);

    fid = fopen(json_filename, 'w');
    if fid == -1, error('Cannot create json file'); end
    fwrite(fid, json_string, 'char');
    fclose(fid);

    disp(['succesfully dumped gnsf model into ', json_filename])
end


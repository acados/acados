function [ output ] = proto_gnsf(s, opts, x0, u0)
% nlf_proto is a prototype to simulate a system with NLF structure + DAE 
% and a Linear Output system as proposed sec.3.3 in Jonathan Freys Master
% thesis draft
% Author: Jonathan Frey
%% Initialization
if isfield(opts, 'n_steps')
    n_steps = opts.n_steps;
else
    n_steps = 1;
end
if isfield(opts, 'S_forw_in')
    S_forw = S_forw_in;
    opts.forw_sens_mode = 1;
else
    S_forw = eye(s.nx, s.nx + s.nu);
end
if isfield(opts, 'adj_seed')
    lambda_w = opts.adj_seed;
else
    lambda_w = eye( s.nx + s.nu,1);
end
time_forw = 0;
time_adj  = 0;
% eps = 1E-11;
q = length(opts.b_butcher);
nff = q * s.n_out;

Z_traj  = zeros( s.nz * q , n_steps);
K1_traj = zeros( s.nx1 * q, n_steps);
x1_traj = zeros( s.nx1 * q, n_steps);
x0_traj = zeros( s.nx * (n_steps+1),1);
x0_traj(1:s.nx,1) = x0;
fftraj  = zeros( s.n_out * q, n_steps);
yytraj  = zeros( s.ny  * q, n_steps);
f_LO_traj = zeros( s.nx2 * q * n_steps, 1 + 2*s.nx1 + s.nu + s.nz);
J_G2_K1 = zeros(s.nx2 * q, s.nx1 * q);
Phi_inc_dy_uhat = zeros(s.n_out * q, 1 + s.ny + s.nuhat);

aux_G2_ff = zeros( s.nx2 * q, q * s.n_out); % J_G2_Z * s.ZZf;
aux_G2_x1 = zeros( s.nx2 * q, s.nx1); % J_G2_Z * s.ZZx
aux_G2_u  = zeros( s.nx2 * q, s.nu ); % J_G2_Z * s.ZZu

r_J_r_ffx1u = zeros(nff, nff + s.nx1 + s.nu +1);
yyu = s.YYu * u0;
uhat = s.L_u * u0;
for ss = 1:n_steps
    % Initialization inside
    x0_1 = x0_traj( (1:s.nx1)      + s.nx *(ss-1) );
    x0_2 = x0_traj( (1+s.nx1:s.nx) + s.nx *(ss-1) );
    yyss = s.YYx * x0_1 + yyu;
    %% Simulation of NLF & DAE system - Newton scheme on ff, Z
    for iter = 1:opts.max_newton    % while max(abs(res_val)) > eps %&& i_newton < opts.max_newton % || i_newton < 10
        % update yy
        yytraj(:, ss) =  yyss + s.YYf * fftraj(:,ss);
        % eval residual, dr_dyy
        r_J_r_ffx1u(:,2:1+nff) = eye(nff); % dr_dff
        for ii = 1:q % eval Phi, Phi_dy;
            ind_y = index(s.ny, ii);
            ind_Phi = index(s.n_out, ii);
            Phi_inc_dy_uhat(ind_Phi,1:(1+s.ny)) = full(s.phi_fun_jac_y(yytraj(ind_y,ss),uhat));
            r_J_r_ffx1u(ind_Phi, 2:nff+1) = r_J_r_ffx1u(ind_Phi, 2:nff+1) - (Phi_inc_dy_uhat(ind_Phi,2:1+s.ny) * s.YYf(ind_y, :)); %build J_r_ff
        end
        r_J_r_ffx1u(:,1) = fftraj(:,ss) - Phi_inc_dy_uhat(:,1);  % compute resval
        delta_ff =  - r_J_r_ffx1u(:,2: nff+1) \ r_J_r_ffx1u(:,1); % solve Newton linear equation system
        fftraj(:,ss) = fftraj(:,ss) + delta_ff;
    end % end Newton iteration
    K1_traj(:,ss) = s.KKf * fftraj(:,ss) + s.KKu * u0 + s.KKx * x0_1;
    Z_traj(:,ss) = s.ZZf * fftraj(:,ss) + s.ZZu * u0 + s.ZZx * x0_1;
    %% build x1_stage_val
    x1_traj(:,ss) = repmat( x0_1, q, 1); %x1_stage_val = repmat( x0_1, q, 1);
    for ii = 1:q %q
        ind_x1_ii = index(s.nx1, ii);
        for jj = 1:q
            ind_x1_jj = index(s.nx1, jj);
            x1_traj(ind_x1_ii,ss) = x1_traj(ind_x1_ii,ss) +...
                opts.A_dt(ii,jj) * K1_traj(ind_x1_jj,ss);
        end
    end
    %% Simulation of Linear Output System:
    % Build G2_value in 1st column of f_LO_traj and get f_LO derivatives in the subsequent columns
    ind_f_LO = index( s.nx2 * q, ss);
    f_LO_traj( ind_f_LO, 1 ) = repmat( - s.ALO * x0_2, q, 1);
    for ii = 1:q
        ind_f_LO_ii = index( s.nx2, ii) + q * s.nx2 * (ss-1);
        ind_x1 = index( s.nx1, ii);
        ind_z  = index( s.nz,  ii);
        f_LO_traj( ind_f_LO_ii, :) = f_LO_traj( ind_f_LO_ii, :) - ...
            full(s.f_lo_fun_jac_x1k1uz(x1_traj(ind_x1, ss), K1_traj(ind_x1,ss), Z_traj(ind_z,ss), u0));
    end
    K2_val = - s.M2inv * f_LO_traj(ind_f_LO,1);
    %% Get simulation result
    x0_traj((1:s.nx) + ss * s.nx) = x0_traj((1:s.nx) + (ss-1) * s.nx); % copy last value
    for ii = 1:q
        ind_x1 = index(s.nx1, ii);
        ind_x2 = index(s.nx2, ii);
        x0_traj((1:s.nx) + ss * s.nx) = x0_traj((1:s.nx) + ss * s.nx) + opts.b_dt(ii) * [K1_traj( ind_x1, ss); K2_val( ind_x2)];
    end
    %% Direct sensitivity propagation
    % eval Jacobian mem.J_PHI_yy, J_PHI_u
    if opts.forw_sens_mode
        tic
        %% eval J_r_ff & J_r_x1u (rhs for sensitivity propagation)
        r_J_r_ffx1u(:,2:1+nff) = eye(nff); % dr_dff
        ind_jac = 2:s.ny+s.nuhat+1;
        ind_jac_y = 2:s.ny+1;
        for ii = 1:q % eval Phi, Phi_dy;
            ind_y = index(s.ny, ii);
            ind_Phi = index(s.n_out, ii);
            Phi_inc_dy_uhat(ind_Phi,ind_jac) = full(s.phi_jac_y_uhat(yytraj(ind_y,ss), uhat));
%             keyboard
            r_J_r_ffx1u(ind_Phi, 2:nff+1) = r_J_r_ffx1u(ind_Phi, 2:nff+1) - Phi_inc_dy_uhat(ind_Phi,ind_jac_y) * s.YYf(ind_y, :); %build J_r_ff
            r_J_r_ffx1u(ind_Phi, 2+nff:1+nff+s.nx1) = - Phi_inc_dy_uhat(ind_Phi,ind_jac_y) * s.YYx(ind_y, :); % build J_r_x1
            r_J_r_ffx1u(ind_Phi, 2+nff+s.nx1:1+nff+s.nx1+s.nu) = - Phi_inc_dy_uhat(ind_Phi,ind_jac_y) * s.YYu(ind_y, :) - Phi_inc_dy_uhat(ind_Phi, s.ny+2:s.ny+s.nuhat+1) * s.L_u; % build J_r_u
        end
        ind_ff = index( s.n_out * q, 1);
        dff_dx1u  = - r_J_r_ffx1u( ind_ff, ind_ff+1) \ r_J_r_ffx1u( :, s.n_out * q + 2: s.n_out*q+ s.nx1 + s.nu+1);
        
        dff_dx1u
        %% expand sensitivities
        %dK1_dwn = [s.KKx + s.KKf * dffZ_dx1u( 1:s.n_out*n_stages, 1:s.nx1), zeros(s.nx1 * n_stages, s.nx2),...
        %             s.KKu + s.KKf * dffZ_dx1u( 1:s.n_out*n_stages, 1+s.nx1:s.nx1+s.nu)];
        % Sensitivities LOS  f_LO_inc_J_x1uz_fun
        dK1_dx1 = s.KKx + s.KKf * dff_dx1u( :, 1:s.nx1);
        dK1_du  = s.KKu + s.KKf * dff_dx1u( :, 1+s.nx1:s.nx1+s.nu);
        dZ_dx1  = s.ZZx + s.ZZf * dff_dx1u(: , 1: s.nx1);
        dZ_du   = s.ZZu + s.ZZf * dff_dx1u(: , 1+s.nx1:s.nx1+s.nu);
        % build J_G2_wn, J_G2_K1
        for ii = 1: q
            ind_x2 = index(s.nx2, ii);
            ind_f_LO_ii = index(s.nx2, ii) + (ss-1) * q * s.nx2;
            ind_z  = index(s.nz , ii);
            ind_x1_ii = index( s.nx1, ii);
            for jj = 1:q
               ind_x1_jj = index( s.nx1, jj); 
               J_G2_K1(  ind_x2, ind_x1_jj) = f_LO_traj( ind_f_LO_ii , 2:s.nx1+1) * opts.A_dt(ii,jj);
            end
            J_G2_K1( ind_x2, ind_x1_ii) = J_G2_K1( ind_x2, ind_x1_ii) + f_LO_traj(ind_f_LO_ii, s.nx1+2: 2*s.nx1+1);
            aux_G2_x1( ind_x2, :) = f_LO_traj( ind_f_LO_ii, 2+2*s.nx1: 2*s.nx1+s.nz +1) * dZ_dx1(ind_z, :);
            aux_G2_u( ind_x2, :) =  f_LO_traj( ind_f_LO_ii, 2+2*s.nx1: 2*s.nx1+s.nz +1) * dZ_du(ind_z, :);
        end
        dK2_dx1(:,1:s.nx1) = -s.M2inv * ( f_LO_traj(ind_f_LO, 2:s.nx1+1) + ...
            J_G2_K1 * dK1_dx1 +...
            aux_G2_x1); % dK2_dx1
        dK2_du = -s.M2inv * ( f_LO_traj(ind_f_LO, 2+2*s.nx1+s.nz: 2*s.nx1+ s.nu+s.nz+1) + J_G2_K1 * dK1_du +...
            aux_G2_u); %dK2_du
        dxf_dwn = eye( s.nx, s.nx + s.nu); % d_psi_wn
        for ii=1:q
            ind_x1 = index(s.nx1, ii);
            ind_x2 = index(s.nx2, ii);
            dxf_dwn(:, 1:s.nx1) = dxf_dwn(:, 1:s.nx1) + opts.b_dt(ii) * ...
                    [dK1_dx1(ind_x1,:); dK2_dx1(ind_x2,:)];  %dxf_dx1
            dxf_dwn(1+s.nx1:s.nx, s.nx1+1:s.nx) = dxf_dwn(1+s.nx1:s.nx, s.nx1+1:s.nx) + ...
                    opts.b_dt(ii) * s.dK2_dx2(ind_x2,:); %dxf_dx2
            dxf_dwn(:, s.nx+1:s.nx+s.nu) = dxf_dwn(:, s.nx+1:s.nx+s.nu) + ...
                opts.b_dt(ii) * [dK1_du(ind_x1,:); dK2_du(ind_x2,:)]; %dxf_du
        end
        dxf_dw0 = [dxf_dwn(1:s.nx, 1:s.nx) * S_forw(1:s.nx, 1:s.nx), ...
            dxf_dwn(1:s.nx, 1:s.nx) * S_forw(:, s.nx+1: s.nx+s.nu) + dxf_dwn(:, s.nx+1: s.nx+ s.nu)]; % multiply with forward seed
        output.forw_sensi = dxf_dw0;
        % update for next step:
        S_forw = output.forw_sensi;
        time_forw = time_forw + toc;
    end
end
output.xf = x0_traj( (1:s.nx) + n_steps * s.nx); % simulation result
for i = 1:s.nz
    for ii =1:opts.n_stages
        Z_i(ii) = Z_traj(s.nz*(ii-1)+i,1);
    end
    output.Zf(i) = neville(0, q-1, opts.c_butcher, Z_i);
end
% z_poly = polyfit(opts.c_butcher, Ztraj(:,1),3);
if opts.adj_sens_mode
    tic
    for ss = n_steps:-1:1
    %% adjoint sensitivity propagation
    x0_1 = x0_traj( (1:s.nx1)      + s.nx *(ss-1) );
    % build J_G2_wn, J_G2_K1
    ind_f_LO = index( s.nx2 * q, ss);
    for ii = 1: q
        ind_x2 = index(s.nx2, ii);
        ind_z  = index(s.nz , ii);
        ind_f_LO_ii = index(s.nx2, ii) + (ss-1) * q * s.nx2;
        aux_G2_ff( ind_x2, : ) = f_LO_traj( ind_f_LO_ii, 2+2*s.nx1: 2*s.nx1+s.nz +1) * s.ZZf(ind_z,:); 
        aux_G2_x1( ind_x2, : ) = f_LO_traj( ind_f_LO_ii, 2+2*s.nx1: 2*s.nx1+s.nz +1) * s.ZZx(ind_z,:);
        aux_G2_u ( ind_x2, : ) = f_LO_traj( ind_f_LO_ii, 2+2*s.nx1: 2*s.nx1+s.nz +1) * s.ZZu(ind_z,:);
        for jj = 1:q
           ind_x1_jj = index( s.nx1, jj); 
           J_G2_K1(  ind_x2, ind_x1_jj) = f_LO_traj( ind_f_LO_ii , 2:s.nx1+1) * opts.A_dt(ii,jj);
        end
        J_G2_K1( ind_x2, ind_x1_ii ) = J_G2_K1( ind_x2, ind_x1_ii ) + f_LO_traj(ind_f_LO_ii, 2+s.nx1:1+2*s.nx1);
    end
    dK2_dff = - s.M2inv * (J_G2_K1 * s.KKf + aux_G2_ff); % J_G2_Z * s.ZZf);
    dK2_dx1 = - s.M2inv * (J_G2_K1 * s.KKx + aux_G2_x1 + f_LO_traj( ind_f_LO, 2:s.nx1+1));
    dK2_du  = - s.M2inv * (f_LO_traj( ind_f_LO, 2+2*s.nx1+s.nz:2*s.nx1+s.nu+s.nz+1)  + J_G2_K1 * s.KKu + aux_G2_u);
    
    dPsi_dff = zeros(s.nx, s.n_out * q);  % initialize in each stage, because we just add stuff
    dPsi_dx  = eye( s.nx, s.nx);
    dPsi_du  = zeros( s.nx, s.nu );
    ind_Psi1 = 1:s.nx1;
    ind_Psi2 = s.nx1+1 : s.nx;
    for ii = 1:q
        ind_x1 = index(s.nx1,ii);
        ind_x2 = index(s.nx2, ii);
        dPsi_dff(ind_Psi1,:) = dPsi_dff(ind_Psi1,:) + opts.b_dt(ii) * s.KKf( ind_x1,:); %Psi1
        dPsi_dx( ind_Psi1, ind_Psi1)  = dPsi_dx(ind_Psi1, ind_Psi1) + opts.b_dt(ii) * s.KKx(ind_x1,:); % dPsi1_dx1
        dPsi_du( ind_Psi1, :)  = dPsi_du(1:s.nx1, :) + opts.b_dt(ii) * s.KKu(ind_x1,:); %dPsi1_du
        dPsi_dff(ind_Psi2, :) = dPsi_dff(ind_Psi2,:) + opts.b_dt(ii) * dK2_dff(ind_x2,:); %Psi2
        dPsi_dx(ind_Psi2, :) = dPsi_dx(ind_Psi2, :) + opts.b_dt(ii) * [ dK2_dx1(ind_x2,:), s.dK2_dx2(ind_x2,:)];
        dPsi_du( ind_Psi2, :) = dPsi_du( ind_Psi2,:) + opts.b_dt(ii) * dK2_du(ind_x2,:); %dPsi2_du
    end
%     r_J_r_ffx1u(:, 2 : s.n_out*q + s.nx1 + s.nu+1) = full(s.jac_res_ffx1u_fun(fftraj(:,ss), x0_1, u0));

    %% eval J_r_ff & J_r_x1u (rhs for sensitivity propagation
    r_J_r_ffx1u(:,2:1+nff) = eye(nff); % dr_dff
    ind_jac = 2:s.ny+s.nuhat+1;
    ind_jac_y = 2:s.ny+1;
    for ii = 1:q % eval Phi, Phi_dy;
        ind_y = index(s.ny, ii);
        ind_Phi = index(s.n_out, ii);
        Phi_inc_dy_uhat(ind_Phi,ind_jac) = full(s.phi_jac_y_uhat(yytraj(ind_y,ss), uhat));

        r_J_r_ffx1u(ind_Phi, 2:nff+1) = r_J_r_ffx1u(ind_Phi, 2:nff+1) - Phi_inc_dy_uhat(ind_Phi,ind_jac_y) * s.YYf(ind_y, :); %build J_r_ff
        r_J_r_ffx1u(ind_Phi, 2+nff:1+nff+s.nx1) = - Phi_inc_dy_uhat(ind_Phi,ind_jac_y) * s.YYx(ind_y, :); % build J_r_x1
%         r_J_r_ffx1u(ind_Phi, 2+nff+s.nx1:1+nff+s.nx1+s.nu) = - Phi_inc_dy_uhat(ind_Phi,ind_jac_y) * s.YYu(ind_y, :); % build J_r_u
        r_J_r_ffx1u(ind_Phi, 2+nff+s.nx1:1+nff+s.nx1+s.nu) = - Phi_inc_dy_uhat(ind_Phi,ind_jac_y) * s.YYu(ind_y, :) - Phi_inc_dy_uhat(ind_Phi, s.ny+2:s.ny+s.nuhat+1) * s.L_u; % build J_r_u
    end
    ind_ff = index( s.n_out * q, 1);
    lambda_ff = - (transpose(r_J_r_ffx1u(ind_ff,ind_ff+1))) \ ( transpose(dPsi_dff)  * lambda_w(1:s.nx,:));
    lambda_w = transpose([dPsi_dx, dPsi_du; zeros(s.nu,s.nx), eye(s.nu,s.nu)])  * lambda_w;
    lambda_w = lambda_w + transpose( ...
         [ r_J_r_ffx1u(:, s.n_out*q+2 : q*s.n_out+s.nx1+1), zeros( s.n_out*q, s.nx2), r_J_r_ffx1u(:, s.n_out*q+s.nx1+2 : q*s.n_out + s.nx1 + s.nu+1)]) * lambda_ff;
    end
    time_adj = toc;
end
output.adj_sensi = lambda_w;
output.time_forw = time_forw;
output.time_adj  = time_adj;
% keyboard
end
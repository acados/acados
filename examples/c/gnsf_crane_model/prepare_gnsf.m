function [ s ] = prepare_gnsf2(s,opts)
tic
    import casadi.*
    nx1 = s.nx1; nx = s.nx; nu = s.nu; nz = s.nz; n_out = s.n_out;
    I_stages = eye( opts.n_stages );
    %generate submatrices
    E11 = s.E(1:nx1, 1:nx1);
    E12 = s.E(1:nx1, 1+nx1:nx1+nz);
    E21 = s.E(1+nx1:nx1+nz, 1:nx1);
    E22 = s.E(1+nx1:nx1+nz, 1+nx1:nx1+nz);
    
    A1 = s.A(1:nx1, :); A2 = s.A(nx1+1:nx1+nz, :);
    B1 = s.B(1:nx1, :); B2 = s.B(nx1+1:nx1+nz, :);
    C1 = s.C(1:nx1, :); C2 = s.C(nx1+1:nx1+nz, :);
    
    % generate fat matrices
    AA1 = repmat(A1,opts.n_stages,1);
    AA2 = repmat(A2,opts.n_stages,1);
    BB1 = repmat(B1,opts.n_stages,1);
    BB2 = repmat(B2,opts.n_stages,1);
    CC1 = kron(I_stages,C1);
    CC2 = kron(I_stages,C2);
    DD1 = -kron(I_stages, E12);
    DD2 = opts.dt * kron(opts.A_butcher, A2) - kron(I_stages, E21);
    EE1 = kron(I_stages, E11) - opts.dt * kron(opts.A_butcher, A1);
    EE2 = kron(I_stages, E22);
    
    EE1inv = inv(EE1);
    EE2inv = inv(EE2);
    
    if 0%nx1 < nz
        % needed to get K-values from f,u,x
        PP1 = inv(eye(nx1*opts.n_stages) - (EE1 \ DD1) * EE2inv  * DD2); 
        PP2 = (EE1 \ DD1) * EE2inv;
        KKf = PP1 * ( PP2 * CC2 + EE1 \ CC1);
        KKu = PP1 * ( PP2 * BB2 + EE1 \ BB1);
        KKx = PP1 * ( PP2 * AA2 + EE1 \ AA1);

        % to get Z-values
        ZZf = EE2 \ (DD2 * KKf + CC2);
        ZZu = EE2 \ (DD2 * KKu + BB2);
        ZZx = EE2 \ (DD2 * KKx + AA2);
    else
        QQ1 = inv(eye(nz*opts.n_stages) - (EE2 \ DD2) * (EE1 \DD1));
        QQ2 = EE2 \ DD2 * EE1inv;
        
        ZZf = QQ1 * ( QQ2 * CC1 + EE2 \ CC2);
        ZZu = QQ1 * ( QQ2 * BB1 + EE2 \ BB2);
        ZZx = QQ1 * ( QQ2 * AA1 + EE2 \ AA2);
        
        % to get K-values
        KKf = EE1 \ (DD1 * ZZf + CC1);
        KKu = EE1 \ (DD1 * ZZu + BB1);
        KKx = EE1 \ (DD1 * ZZx + AA1);
                
    end

    % to get yy-values:
    LLZ = kron(I_stages, s.L_z);
    LLx = repmat(s.L_x, opts.n_stages, 1);
    LLK = kron(opts.dt * opts.A_butcher, s.L_x) + kron(I_stages, s.L_xdot);
%     LLu = repmat(s.L_u, opts.n_stages, 1);

    YYx = LLK * KKx + LLZ * ZZx + LLx;
    YYu = LLK * KKu + LLZ * ZZu;
    YYf = LLK * KKf + LLZ * ZZf;
%     keyboard
    s.YYx = YYx;
    s.YYu = YYu;
    s.YYf = YYf;
      
    % needed for linear output system:
    M2 = eye( s.nx2 * opts.n_stages);
    M2 =  M2 - kron( opts.dt * opts.A_butcher, s.ALO);
    M2inv = M2^-1;

    % put precomputed matrices in struct
    s.KKf = KKf;
    s.KKu = KKu;
    s.KKx = KKx;
    
    s.ZZf = ZZf;
    s.ZZu = ZZu;
    s.ZZx = ZZx;
    
    s.M2inv = M2inv;
    
    s.dK2_dx2 = M2inv * repmat(s.ALO,opts.n_stages,1);
    disp(['time in prepare_nlf = ',num2str(toc)]);
%     keyboard
end
function [ H, h, G, g, l, u ] = read_ocp_qp( folder )

if ~ischar(folder)
    error('Input must be string');
end

if ~exist(folder, 'dir')
    error('Input must be valid path');
end

% append trailing `/' if needed
if folder(end) ~= '/'
    folder = [folder, '/'];
end

stage_exists = exist([folder, 'Am0.txt'], 'file');
H = [];
h = [];
G = [];
g = [];
l = [];
u = [];
stage = 0;
while stage_exists
    Q = load([folder, 'Qm', num2str(stage), '.txt']);
    nx = size(Q, 1);
    S = load([folder, 'Sm', num2str(stage), '.txt']);
    nu = size(S, 1);
    R = load([folder, 'Rm', num2str(stage), '.txt']);
    q = load([folder, 'qv', num2str(stage), '.txt']);
    r = load([folder, 'rv', num2str(stage), '.txt']);
    H = blkdiag(H, [Q, S.'; S, R]);
    h = [h; q; r];
    
    A = load([folder, 'Am', num2str(stage), '.txt']);
    B = load([folder, 'Bm', num2str(stage), '.txt']);
    b = load([folder, 'bv', num2str(stage), '.txt']);
    G = [G, zeros(stage*nx, (stage > 0)*(nx+nu)); zeros(nx, stage*(nx+nu)), A, B, -eye(nx)];
    g = [g; b];
    
    lb_stage = load([folder, 'lb', num2str(stage), '.txt']);
    ub_stage = load([folder, 'ub', num2str(stage), '.txt']);
    idxb_stage = floor(load([folder, 'idxb', num2str(stage), '.txt']));

    l_stage = -inf*ones(nx+nu, 1);
    u_stage = +inf*ones(nx+nu, 1);
    
    l_stage(idxb_stage+1) = lb_stage;
    u_stage(idxb_stage+1) = ub_stage;
    
    l = [l; l_stage];
    u = [u; u_stage];
    
    stage = stage+1;
    stage_exists = exist([folder, 'Am', num2str(stage), '.txt'], 'file');
end

Q = load([folder, 'Qm', num2str(stage), '.txt']);
nx = size(Q, 1);
S = load([folder, 'Sm', num2str(stage), '.txt']);
nu = size(S, 1);
R = load([folder, 'Rm', num2str(stage), '.txt']);
q = load([folder, 'qv', num2str(stage), '.txt']);
r = load([folder, 'rv', num2str(stage), '.txt']);
H = blkdiag(H, [Q, S.'; S, R]);
h = [h; q; r];

l = [l; load([folder, 'lb', num2str(stage), '.txt'])];
u = [u; load([folder, 'ub', num2str(stage), '.txt'])];

end


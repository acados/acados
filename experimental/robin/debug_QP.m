clear

file = fopen('QP_data.txt','r');
num_arrays_to_read = 10;
arrays_to_read = {};
for i=1:num_arrays_to_read
    instr=fgets(file);
    temp = sscanf(instr,'%g ');
    arrays_to_read{i} = temp;
end
N = 20;
nx = 8;
nu = 3;
numvars = N*nu;
[H,f,A,lbU,ubU,lbA,ubA,C,d,D] = arrays_to_read{:};
H = reshape(H,numvars,numvars);
A = reshape(A,N*(nx+nx+nu)+nx+nu,numvars);
C = reshape(C,N*(nx),numvars);
D = reshape(D,(N+1)*(nx+nu),numvars);

% test condensing Andersson2013
% Ad = [1,1;0,1];
% Bd = [0.5;1];
% G = kron(eye(N),Bd);
% for i=0:N-1
%     G_block = G(i*nx+1:(i+1)*nx,i*nu+1:(i+1)*nu);
%     for k=i+1:N-1
%         G_block = Ad*G_block;
%         G(k*nx+1:(k+1)*nx,i*nu+1:(i+1)*nu) = G_block;
%     end
% end

u = quadprog(H,f,[A;-A],[ubA;-lbA],[],[],lbU,ubU);

x = C*u+d;
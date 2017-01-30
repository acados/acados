
% generate a toy ocp_qp_in
clc;

d  = 'qp_data';
if ~isdir(d)
    mkdir(d);
end

BOUND  = 5;

N  = 2;
nx = [3 3 3];
nu = [2 2];
nb = nx+[nu 0];
nc = [1 1 1];
x0 = [0.3; 0.3; 0.2];

saveIntVector(N, 'N', d);
saveIntVector(nx, 'nx', d);
saveIntVector(nu, 'nu', d);
saveIntVector(nb, 'nb', d);
saveIntVector(nc, 'nc', d);
saveDoubleVector(x0, 'x0', d);

for ii = 1:N
   A  = ii*ones(nx(ii+1), nx(ii));
   B  = (ii+N)*ones(nx(ii+1), nu(ii));
   b  = (ii+4)*ones(nx(ii+1),1);
   A(1,2) = 5;
   
   Q = rand(nx(ii)) + 10*eye(nx(ii)); Q = 0.5*(Q+Q');
   R = rand(nu(ii)) + 10*eye(nu(ii)); R = 0.5*(R+R');
   S = rand(nu(ii), nx(ii));
   
   q = rand(nx(ii),1);
   r = rand(nu(ii),1);
   
   lb   = -BOUND*ones(nx(ii)+ nu(ii),1);
   ub   = BOUND*ones(nx(ii)+ nu(ii),1);
   idxb = (1:nx(ii)+ nu(ii))-1;

   lc = -2*ones(nc(ii),1);
   uc = 2*ones(nc(ii),1);
   Cx = ones(nc(ii), nx(ii));
   Cu = ones(nc(ii), nu(ii));
   
   saveDoubleMatrix(A, ['A' num2str(ii-1)], d);
   saveDoubleMatrix(B, ['B' num2str(ii-1)], d);
   saveDoubleVector(b, ['b' num2str(ii-1)], d);

   saveDoubleMatrix(Q, ['Q' num2str(ii-1)], d);
   saveDoubleMatrix(R, ['R' num2str(ii-1)], d);
   saveDoubleMatrix(S, ['S' num2str(ii-1)], d);

   saveDoubleVector(q, ['q' num2str(ii-1)], d);
   saveDoubleVector(r, ['r' num2str(ii-1)], d);
   
   saveDoubleVector(lb, ['lb' num2str(ii-1)], d);
   saveDoubleVector(ub, ['ub' num2str(ii-1)], d);
   saveIntVector(idxb, ['idxb' num2str(ii-1)], d);
   
   saveDoubleMatrix(Cx, ['Cx' num2str(ii-1)], d);
   saveDoubleMatrix(Cu, ['Cu' num2str(ii-1)], d);
   saveDoubleVector(lc, ['lc' num2str(ii-1)], d);
   saveDoubleVector(uc, ['uc' num2str(ii-1)], d);
end
lb   = -0.45*ones(nx(ii),1);
ub   = 0.68*ones(nx(ii),1);
idxb = (1:nx(ii))-1;
Q = rand(nx(ii)) + 10*eye(nx(ii)); Q = 0.5*(Q+Q');
q = rand(nx(ii),1);
saveDoubleMatrix(Q, ['Q' num2str(N)], d); 
saveDoubleVector(q, ['q' num2str(N)], d);
saveDoubleVector(lb, ['lb' num2str(N)], d);
saveDoubleVector(ub, ['ub' num2str(N)], d);
saveIntVector(idxb, ['idxb' num2str(N)], d);

lc = -0.28;
uc = 0.28;
Cx = [0 1 0];
saveDoubleMatrix(Cx, ['Cx' num2str(N)], d);
saveDoubleVector(lc, ['lc' num2str(N)], d);
saveDoubleVector(uc, ['uc' num2str(N)], d);

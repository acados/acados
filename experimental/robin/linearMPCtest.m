clear

N = 20;
nx = 8;
nu = 3;

A = [0.76272   0.11488   0.00248   0.00002   0.45961   0.01981   0.00025   0.00000 ;
  0.11488   0.76520   0.11490   0.00248   0.01981   0.45987   0.01981   0.00025 ;
  0.00248   0.11490   0.76520   0.11488   0.00025   0.01981   0.45987   0.01981 ;
  0.00002   0.00248   0.11488   0.76272   0.00000   0.00025   0.01981   0.45961 ;
 -0.89941   0.42024   0.01931   0.00025   0.76272   0.11488   0.00248   0.00002 ;
  0.42024  -0.88010   0.42049   0.01931   0.11488   0.76520   0.11490   0.00248 ;
  0.01931   0.42049  -0.88010   0.42024   0.00248   0.11490   0.76520   0.11488 ;
  0.00025   0.01931   0.42024  -0.89941   0.00002   0.00248   0.11488   0.76272 ];

B = [0.11990   0.00252   0.00002 ;
  0.00252   0.11992   0.00252 ;
  0.00002   0.00252   0.11992 ;
  0.00000   0.00002   0.00252 ;
  0.45961   0.01981   0.00025 ;
  0.01981   0.45987   0.01981 ;
  0.00025   0.01981   0.45987 ;
  0.00000   0.00025   0.01981 ];

b = [0.10000 ;
  0.10000 ;
  0.10000 ;
  0.10000 ;
  0.10000 ;
  0.10000 ;
  0.10000 ;
  0.10000];

x0 = [  2.50000 ;
  2.50000 ;
  0.00000 ;
  0.00000 ;
  0.00000 ;
  0.00000 ;
  0.00000 ;
  0.00000 ];

% Cost function
Q = eye(nx);
R = 2*eye(nu);
q = 0.1*ones(nx,1);
r = 0.2*ones(nu,1);
H = blkdiag(kron(eye(N),blkdiag(Q,R)),Q);
f = [repmat([q;r],N,1);q];

C = [-eye(nx),zeros(nx,N*(nx+nu))];
c = [x0;repmat(b,N,1)];
for i=0:N-1
    C = [C;zeros(nx,i*(nx+nu)),A,B,-eye(nx),zeros(nx,(N-i-1)*(nx+nu))];
end


lb = [repmat([-4*ones(nx,1);-0.5*ones(nu,1)],N,1);-4*ones(nx,1)];
lb(end-nx+1:end) = 0;
ub = -lb;


sparse.w = quadprog(H,f,[],[],C,-c,lb,ub);

XU = reshape([sparse.w;zeros(nu,1)],nx+nu,N+1);
sparse.x = vec(XU(1:nx,:));
sparse.u = vec(XU(nx+1:end,1:end-1));

figure(1);clf;plot(XU(1:2:7,:).')

%% Condensing

Abar = [-eye(nx),zeros(nx,(N-1)*nx)];
Bbar = kron(eye(N),B);
for i=1:N-1
    Abar = [Abar;zeros(nx,(i-1)*(nx)),A,-eye(nx),zeros(nx,(N-i-1)*nx)];
end
Qbar = kron(eye(N+1),Q);
Rbar = kron(eye(N),R);
c = [A*x0+b;repmat(b,N-1,1)];
Dx = blkdiag(zeros(N*nx,N*nx),eye(nx));

G = [zeros(nx,N*nu);-Abar\Bbar];
Ge = [eye(nx);-Abar\[A;zeros((N-1)*nx,nx)]];
g = [zeros(nx,1);-Abar\c];

Hbar = Rbar + G.'*Qbar*G;
hbar = repmat(r,N,1) + G.'*(repmat(q,N+1,1)+Qbar*g);
%% Read from file
file = fopen('QP_data.txt','r');
num_arrays_to_read = 10;
arrays_to_read = {};
for i=1:num_arrays_to_read
    instr=fgets(file);
    temp = sscanf(instr,'%g ');
    arrays_to_read{i} = temp;
end
numvars = N*nu;
[H,f,A,lbU,ubU,lbA,ubA,C,d,D] = arrays_to_read{:};
H = reshape(H,numvars,numvars);
A = reshape(A,N*(nx+nx+nu)+nx+nu,numvars);
C = reshape(C,N*(nx),numvars);
D = reshape(D,(N+1)*(nx+nu),numvars);

bigA = [];
for i=0:N-1
    bigA = [bigA;C(i*nx+1:(i+1)*nx,:);D(i*(nx+nu)+1:(i+1)*(nx+nu),:)];
end
bigA = [bigA;D(end-(nx+nu)+1:end,:)];

condensing.u = quadprog(H,hbar,[A;-A],[ubA;-lbA],[],[],lbU,ubU);
condensing.x = C*condensing.u+d;

norm(condensing.u-sparse.u)
norm(condensing.x-sparse.x(nx+1:end))

figure(1);hold on;
condensing.XU = reshape([x0;condensing.x],nx,N+1);
plot(condensing.XU(1:2:7,:).','o')
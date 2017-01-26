clc;
clear all;
close all;

Ts = 0.01;
EXPORT = 1;

DifferentialState x1 x2;
Control u;

xmin = -2; xmax = 2;
umin = -10; umax = 10;

sqp_iter = 50;
%% Differential Equation
    
% Model equations
x = [x1 x2];
f_expl = ode(1,x,u);

%% SIMexport
acadoSet('problemname', 'sim');

numSteps = 100;
sim = acado.SIMexport( Ts );
sim.setModel(f_expl);
sim.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2' );
sim.set( 'NUM_INTEGRATOR_STEPS',        numSteps        );

if EXPORT
    sim.exportCode( 'export_SIM' );
    
    cd export_SIM
    make_acado_integrator('../integrate_system')
    cd ..
end

%% MPCexport
acadoSet('problemname', 'mpc');

N = 13;
ocp = acado.OCP( 0.0, N*Ts, N );

h = [x1 x2 u];
hN = [x1 x2];

W = acado.BMatrix(eye(length(h)));
WN = acado.BMatrix(eye(length(hN)));
 
ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );



ocp.subjectTo( xmin <= x1 <= xmax );
ocp.subjectTo( xmin <= x2 <= xmax );
ocp.subjectTo( umin <= u <= umax );
% ocp.subjectTo( 'AT_END', [x1 theta v1 dtheta] == 0 );

ocp.setModel(f_expl);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
% mpc.set( 'LEVENBERG_MARQUARDT',          1e-5               );
mpc.set( 'INTEGRATOR_TYPE',             'INT_RK4'           );
mpc.set( 'NUM_INTEGRATOR_STEPS',        2*N                 );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);

if EXPORT
    mpc.exportCode( 'export_MPC' );
    
    global ACADO_;
    copyfile([ACADO_.pwd '/../../external_packages/qpoases'], 'export_MPC/qpoases')
    
    cd export_MPC
    make_acado_solver('../acado_MPCstep')
    cd ..
end

%% PARAMETERS SIMULATION
X0 = [0.5 0];
Xref = [0 0];
input.x = zeros(N+1,2);
input.x0 = X0;

Uref = zeros(N,1);
input.u = Uref;

input.y = [repmat(Xref,N,1) Uref];
input.yN = Xref;

input.W = diag([1 1 0.05]);
input.WN = diag([1 1]);

input.shifting.strategy =0;

KKT_MPC = []; INFO_MPC = [];
controls_MPC = [];
state_sim = X0;


for iter = 1:sqp_iter
    tic
    % Solve NMPC OCP
%     input.x0 = state_sim(end,:);
    output = acado_MPCstep(input);
    
    % Save the MPC step
    INFO_MPC = [INFO_MPC; output.info];
    KKT_MPC = [KKT_MPC; output.info.kktValue];
    controls_MPC = [controls_MPC; output.u(1,:)];

    input.x = output.x;
    input.u = output.u;
        
    disp(['current iteration: ' num2str(iter) '   ' char(9) ' (RTI step -- QP error: ' num2str(output.info.status) ',' ' ' char(2) ' KKT val: ' num2str(output.info.kktValue,'%1.2e') ',' ' ' char(2) ' CPU time: ' num2str(round(output.info.cpuTime*1e6)) ' µs)'])
    
end

figure()
subplot(211)
plot(Ts*[0:N],output.x)
xlabel('time [s]')
ylabel('states')
legend('x_1','x_2')
grid on

subplot(212)
stairs(Ts*[0:N-1],output.u)
xlabel('time [s]')
ylabel('u')
grid on

[output.x [output.u;0]]
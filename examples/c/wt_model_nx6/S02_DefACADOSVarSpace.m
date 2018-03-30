%% States
import casadi.*

% States 
x = MX.sym('x',6); 

% Generator Angular Velocity ( rad/s ) 
GEN_agvelSt = x(1);
% Drivetrain torsional Angular Velocity ( rad/s ) 
DT_agvelTorsSt = x(2);
% Generator azimuth angle ( rad ) 
GEN_agSt = x(3);
% Drivetrain torsional angle ( rad ) 
DT_agTorsSt = x(4);
% Blade pitch angle ( rad ) 
BLD_agPtchActSt = x(5);
% Generator torque ( 10kNm ) 
GEN_trqActSt = x(6);

%% Differential State Variables 

% Differntial States 
dx = MX.sym('dx',6); 

% Drivetrain angular acceleration ( rad/s^2 ) 
dotDT_agaccDyn = dx(1);
% Drivetrain torsional angular acceleration ( rad/s^2 ) 
dotDT_agaccDynTors = dx(2);
% Drivetrain angular velocity ( rad/s ) 
dotDT_Dyn_AngVel = dx(3);
% Drivetrain torsional angular velocity ( rad/s ) 
dotDT_Dyn_AngVel_Tors = dx(4);
% pitch dynamics blade PT-1 ( rad/s ) 
dotBLD_Dyn_PtchAct = dx(5);
% Generator torque PT-1 ( Nm/s ) 
dotGEN_Dyn_TrqAct = dx(6);

%% Control inputs 

% Control inputs 
u = MX.sym('u',2); 

% Desired pitch angle for Blade ( rad ) 
BLD_agPtchDes = u(1);
% Desired generator torque ( 10kNm ) 
GEN_trqDes = u(2);

%% Disturbance 

p = MX.sym('p',1); 

% Far upstream effective wind velocity ( m/s ) 
ENV_velEffWnd = p(1);



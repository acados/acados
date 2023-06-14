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

%% States
import casadi.*

% States 
xy = MX.sym('xy',5);
xnoy = MX.sym('xnoy',1);

% Generator Angular Velocity ( rad/s ) 
GEN_agvelSt = xy(1);
% Drivetrain torsional Angular Velocity ( rad/s ) 
DT_agvelTorsSt = xy(2);
% Generator azimuth angle ( rad ) 
GEN_agSt = xy(3);
% Drivetrain torsional angle ( rad ) 
DT_agTorsSt = xy(4);
% Blade pitch angle ( rad ) 
BLD_agPtchActSt = xy(5);
% Generator torque ( 10kNm ) 
GEN_trqActSt = xnoy(1);

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



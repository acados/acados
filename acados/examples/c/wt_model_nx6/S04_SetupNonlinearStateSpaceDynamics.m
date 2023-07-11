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

%% Set-up nonlinear state-space model of wind turbine

 % Rotor azimuth angle for Blade 1 ( rad ) 
BLD1_agAzi = DT_agTorsSt+GEN_agSt ;

 % Rotor azimuth angle for Blade 2 ( rad ) 
BLD2_agAzi = p_6+DT_agTorsSt+GEN_agSt ;

 % Rotor azimuth angle for Blade 3 ( rad ) 
BLD3_agAzi = p_7+DT_agTorsSt+GEN_agSt ;

 % Blade 1 effective wind speed ( m/s ) 
BLD1_velEff = ((p_3*cos(BLD1_agAzi)+p_4)^(0.2))*(ENV_velEffWnd) ;

 % Blade 2 effective wind speed ( m/s ) 
BLD2_velEff = ((p_3*cos(BLD2_agAzi)+p_4)^(0.2))*(ENV_velEffWnd) ;

 % Blade 3 effective wind speed ( m/s ) 
BLD3_velEff = ((p_3*cos(BLD3_agAzi)+p_4)^(0.2))*(ENV_velEffWnd) ;

 % Rotational speed blade 1 ( - ) 
BLD1_agvel = DT_agvelTorsSt+GEN_agvelSt ;

 % Rotational speed blade 2 ( - ) 
BLD2_agvel = DT_agvelTorsSt+GEN_agvelSt ;

 % Rotational speed blade 3 ( - ) 
BLD3_agvel = DT_agvelTorsSt+GEN_agvelSt ;

 % Tip speed ratio blade 1 ( - ) 
BLD1_velTipRat = p_5*BLD1_agvel/BLD1_velEff ;

 % Tip speed ratio blade 2 ( - ) 
BLD2_velTipRat = p_5*BLD2_agvel/BLD2_velEff ;

 % Tip speed ratio blade 3 ( - ) 
BLD3_velTipRat = p_5*BLD3_agvel/BLD3_velEff ;

 % Individual blade 1 pitch angle ( rad ) 
BLD1_agPtch = BLD_agPtchActSt ;

 % Individual blade 2 pitch angle ( rad ) 
BLD2_agPtch = BLD_agPtchActSt ;

 % Individual blade 3 pitch angle ( rad ) 
BLD3_agPtch = BLD_agPtchActSt ;

 % Tangential aerodynamic force on Blade 1 ( N ) 
BLD1_frTanAero = p_2*BLD1_velEff^p_1*splineCMBL([BLD1_agPtch,BLD1_velTipRat]) ;

 % Tangential aerodynamic force on Blade 2 ( N ) 
BLD2_frTanAero = p_2*BLD2_velEff^p_1*splineCMBL([BLD2_agPtch,BLD2_velTipRat]) ;

 % Tangential aerodynamic force on Blade 3 ( N ) 
BLD3_frTanAero = p_2*BLD3_velEff^p_1*splineCMBL([BLD3_agPtch,BLD3_velTipRat]) ;



%% Explicit Nonlinear State-Space Model
 fe = [ ... 
% Drivetrain angular acceleration ( rad/s^2 ) 
 (DT_agTorsSt*p_14+DT_agvelTorsSt*p_13-GEN_trqActSt*p_12)/(p_10+p_11);
% Drivetrain torsional angular acceleration ( rad/s^2 ) 
 -(BLD1_frTanAero*p_8*p_9+BLD2_frTanAero*p_8*p_9+BLD3_frTanAero*p_8*p_9+DT_agTorsSt*p_14+DT_agvelTorsSt*p_13)/p_10;
% Drivetrain angular velocity ( rad/s ) 
 -GEN_agvelSt*p_8;
% Drivetrain torsional angular velocity ( rad/s ) 
 -DT_agvelTorsSt*p_8;
% pitch dynamics blade PT-1 ( rad/s ) 
 p_15*(-BLD_agPtchActSt+BLD_agPtchDes);
% Generator torque PT-1 ( Nm/s ) 
 p_16*(-GEN_trqActSt+GEN_trqDes);
 ];
 


%% Implicit Nonlinear State-Space Model
fi = dx - fe;

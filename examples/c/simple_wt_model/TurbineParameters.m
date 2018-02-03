function [Parameter] = TurbineParameters


%% General          
Parameter.General.rho                   = 1.225;              % [kg/m^3]  air density

Parameter.Turbine.R                     = 63;

Parameter.Turbine.i                     = 1/97;%1/97;
Parameter.Turbine.J                     = 4.0589e+07;%40.47*10^6;




Parameter.Turbine.k_Te                  = 1.9170e+06; 
Parameter.Turbine.m_Te                  = 4.7218e+05;
Parameter.Turbine.c_Te                  = 1.8834e+04;                   % [kg/s]	tower structual damping (sigma=C_T/(2M_T)=D*w_0, w_0=f_0*2*pi, [Gasch] p.294)

%FILL IN DATA DOWN THERE


%% CPC
Parameter.CPC.Omega_rated               = 1.2671;  % [rad/s]

Parameter.PitchActuator.theta_dot_max   = 5;
Parameter.PitchActuator.theta_max       = 20;  %[deg]
Parameter.PitchActuator.theta_min       =  0;  %[deg]
   

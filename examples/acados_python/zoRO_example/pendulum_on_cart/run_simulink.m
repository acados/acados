% !mkdir -p zoro_test
!cp simulink_model_closed_loop.slx c_generated_code/
clear all;

nx = 4;
nu = 1;
N = 20;
U_max = 40;
ny = nx+nu;
ny_e = nx;
x0 = [0.0, 0.15*pi, 0.0, 0.0];


% in zoro test
cd c_generated_code/
make_sfun_pendulum
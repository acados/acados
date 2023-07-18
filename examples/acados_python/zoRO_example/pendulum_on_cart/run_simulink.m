clear all;

acados_folder = getenv('ACADOS_INSTALL_DIR');

copyfile(fullfile(acados_folder, 'examples', 'acados_matlab_octave', 'getting_started', 'simulink_model_closed_loop.slx'), fullfile(pwd, 'c_generated_code'))

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
clear all;

system('python minimal_example_zoro.py')
acados_folder = getenv('ACADOS_INSTALL_DIR');

copyfile(fullfile(acados_folder, 'examples', 'acados_matlab_octave', 'getting_started', 'simulink_model_closed_loop.slx'), fullfile(pwd, 'c_generated_code'))

nx = 4;
nu = 1;
N = 20;
U_max = 40;
ny = nx+nu;
ny_e = nx;
x0 = [0.0, 0.15*pi, 0.0, 0.0];

cd c_generated_code/
make_sfun_pendulum
make_sfun_sim_pendulum

out_sim = sim('simulink_model_closed_loop', 'SaveOutput', 'on');

u_values = out_sim.logsout{1}.Values.Data;

disp('successfully ran Simulink zoro closed loop simulation');

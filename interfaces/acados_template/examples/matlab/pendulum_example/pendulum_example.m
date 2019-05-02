clear all
close all
clc

import acados_template.*
import casadi.*

% export model 
model = export_ode_model();
    
% create ocp_nlp object
acados_ocp_nlp = acados_ocp_nlp();

% set model_name 
acados_ocp_nlp.model_name = "pendulum_ode";
Tf = 1.0;
nx = 4;
nu = 1;
ny = nx + nu;
nyN = nx;
N = 10;

% set ocp_nlp_dimensions
nlp_dims     = acados_ocp_nlp.dims;
nlp_dims.nx  = nx; 
nlp_dims.ny  = ny; 
nlp_dims.nyN = nyN; 
nlp_dims.nbx = 0;
nlp_dims.nbu = nu; 
nlp_dims.nu  = 1;
nlp_dims.N   = N;

% set weighting matrices
nlp_cost = acados_ocp_nlp.cost;
Q = eye(4);
Q(1,1) = 1e3;
Q(2,2) = 1e-2;
Q(3,3) = 1e3;
Q(4,4) = 1e-2;

R = eye(1);
R(1,1) = 1e-2;

nlp_cost.W = blkdiag(Q, R); 

Vx = zeros(ny, nx);
Vx(1,1) = 1.0;
Vx(2,2) = 1.0;
Vx(3,3) = 1.0;
Vx(4,4) = 1.0;

nlp_cost.Vx = Vx;

Vu = zeros(ny, nu);
Vu(5,1) = 1.0;
nlp_cost.Vu = Vu;

nlp_cost.WN = Q; 

VxN = zeros(nyN, nx);
VxN(1,1) = 1.0;
VxN(2,2) = 1.0;
VxN(3,3) = 1.0;
VxN(4,4) = 1.0;

nlp_cost.VxN = VxN;

nlp_cost.yref  = zeros(ny, 1);
nlp_cost.yrefN = zeros(nyN, 1);

% setting bounds
Fmax = 80.0;
nlp_con = acados_ocp_nlp.constraints;
nlp_con.lbu = [-Fmax];
nlp_con.ubu = [+Fmax];
nlp_con.x0 = [0.0, 0.0, 3.14, 0.0];
nlp_con.idxbu = [1];

% set constants
acados_ocp_nlp.constants.PI  =  3.1415926535897932;

% set QP solver
acados_ocp_nlp.solver_config.qp_solver = "PARTIAL_CONDENSING_HPIPM";
acados_ocp_nlp.solver_config.hessian_approx = "GAUSS_NEWTON";
acados_ocp_nlp.solver_config.integrator_type = "ERK";

% set prediction horizon
acados_ocp_nlp.solver_config.tf = Tf;
acados_ocp_nlp.solver_config.nlp_solver_type = "SQP";

% set header path
acados_ocp_nlp.acados_include_path = "/usr/local/include";
acados_ocp_nlp.acados_lib_path = "/usr/local/lib";

% dump JSON file
json_string = jsonencode(acados_ocp_nlp);
fid = fopen('acados_ocp_nlp.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, json_string, 'char');
fclose(fid);

opts.generate_hess = 1;
generate_c_code_implicit_ode(model, opts);
generate_c_code_explicit_ode(model);

generate_solver('acados_ocp_nlp.json', '/home/andrea/.acados_t/bin/python3')
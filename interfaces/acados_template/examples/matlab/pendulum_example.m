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

function model = export_ode_model()
    import casadi.*
    import acados_template.*
    model_name = 'pendulum_ode';

    % constants
    M = 1;
    m = 0.1;
    g = 9.81;
    l = 0.8;

    % set up states & controls
    x1      = SX.sym('x1');
    theta   = SX.sym('theta');
    v1      = SX.sym('v1');
    dtheta  = SX.sym('dtheta');
    
    x = vertcat(x1, v1, theta, dtheta);

    % controls
    F = SX.sym('F');
    u = vertcat(F);
    
    % xdot
    x1_dot      = SX.sym('x1_dot');
    theta_dot   = SX.sym('theta_dot');
    v1_dot      = SX.sym('v1_dot');
    dtheta_dot  = SX.sym('dtheta_dot');

    xdot = vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot);
    
    % algebraic variables
    z = [];

    % parameters
    p = [];
    
    % dynamics     
    denominator = M + m - m*cos(theta)*cos(theta);
    f_expl = vertcat(v1, (-m*l*sin(theta)*dtheta*dtheta + m*g*cos(theta)*sin(theta)+F)/denominator, dtheta, (-m*l*cos(theta)*sin(theta)*dtheta*dtheta + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator));
    
    f_impl = xdot - f_expl;
   
    model = acados_dae();

    model.f_impl_expr = f_impl;
    model.f_expl_expr = f_expl;
    model.x = x;
    model.xdot = xdot;
    model.u = u;
    model.z = z;
    model.p = p;
    model.name = model_name;
end


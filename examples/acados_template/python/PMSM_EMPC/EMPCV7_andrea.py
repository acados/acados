from acados_template import *
import acados_template as at
import numpy as nmp
from ctypes import *
import matplotlib
import matplotlib.pyplot as plt
import scipy.linalg
import json

CODE_GEN = 1
COMPILE = 1

FORMULATION = 4         # Tracking MPC
                        # ===================
                        # 0 for hexagon (no terminalset) | 5 hexagon + currentslack| 6 PDC Voltage + Hexagon + torque equation | 7 PDC Voltage + currentslack

                        # Economic MPC
                        # ===================
                        # 1 only with Hexagon | 2  PDC Voltage + Hexagon | 3 only with Hexagon  (but pure nonlinear torque, equal to 1) | 4 PDC Voltage + PDC Terminal Set
                        
i_d_ref =  -125         # Setpoints only valid for Formulation 0, 5, 6 and 7
i_q_ref =    10       
w_val = 2000.0          # do not change in this script, sometimes the values below calculate with a fix w_val = 2000 1/S
tau_wal = 10.0 

# constants
L_d = 107e-6
L_q = 150e-6
R_m = 18.15e-3
K_m = 13.8e-3
N_P = 5.0
u_max = 48
i_max = 155.0          # only for simulation

phi = 0.0

x0Start = nmp.array([0, 0])

N = 2      
Ts_sim = 125e-6

WARMSTART_ITERS = 200

INPUT_REG = 1e-2

# # setpoint MPC with hexagon (current slack)
if FORMULATION == 0:        # (works)
    Weight_TUNING = 1
    Weight_E_TUNING = 1
    SLACK_TUNING   =  0
    SLACK_E_TUNING =  0
    wd = 1
    wq = 1
    Ts = 0.000125
    SLACK_TUNINGHessian =  0

# # setpoint MPC with hexagon
# if FORMULATION == 0:
#     Weight_TUNING = 1e2
#     Weight_E_TUNING = 1e2
#     SLACK_TUNING   =  0
#     SLACK_E_TUNING =  0
#     wd = 1
#     wq = 1
#     Ts = 0.000125
#     SLACK_TUNINGHessian =  0

# full EMPC
if FORMULATION == 4:        # (works)
    Weight_TUNING = 1e-1
    Weight_E_TUNING = 1e-1
    SLACK_TUNING   =  1e3
    SLACK_E_TUNING =  1e3
    wd = 1
    wq = 1
    Ts = 0.000250
    SLACK_TUNINGHessian =  0

# # setpoint MPC with hexagon (current slack)
# if FORMULATION == 5:        # (works)
#     Weight_TUNING = 1e2
#     Weight_E_TUNING = 1e3
#     SLACK_TUNING   =  0
#     SLACK_E_TUNING =  0
#     wd = 1
#     wq = 1
#     Ts = 0.000125
#     SLACK_TUNINGHessian =  0

# if FORMULATION == 6:
#     Weight_TUNING = 500
#     Weight_E_TUNING = 500
#     SLACK_TUNING   =  5
#     SLACK_E_TUNING =  500
#     wd = 1
#     wq = 1
#     Ts = 0.000125
#     SLACK_TUNINGHessian =  0

# if FORMULATION == 7:    # (works)
#     Weight_TUNING = 1e2
#     Weight_E_TUNING = 1e2
#     SLACK_TUNING   =  0
#     SLACK_E_TUNING =  0
#     wd = 1
#     wq = 1
#     Ts = 0.000125
#     SLACK_TUNINGHessian =  0

# Model of the PMSM
#====================================================================
def export_ode_model():

    model_name = 'rsm'

    # set up states 
    i_d = SX.sym('i_d')
    i_q = SX.sym('i_q')
    x = vertcat(i_d, i_q)

    # set up controls 
    u_d = SX.sym('u_d')
    u_q = SX.sym('u_q')
    u = vertcat(u_d, u_q)

    # set up xdot 
    i_d_dot = SX.sym('i_d_dot')
    i_q_dot = SX.sym('i_q_dot')
    xdot = vertcat(i_d_dot, i_q_dot)

    # set up parameters
    omega  = SX.sym('omega')        # speed
    dist_d = SX.sym('dist_d')       # d disturbance
    dist_q = SX.sym('dist_q')       # q disturbance
    tau_des = SX.sym('tau_des')     # Desired tau
    p = vertcat(omega, dist_d, dist_q, tau_des)

    # dynamics     
    fexpl = vertcat(-(R_m/L_d)*i_d + (L_q/L_d)*omega*i_q + u_d/L_d, \
                    -(L_d/L_q)*omega*i_d - (R_m/L_q)*i_q + u_q/L_q - (omega*K_m)/L_q)

    fimpl = vertcat(-i_d_dot-(R_m/L_d)*i_d + (L_q/L_d)*omega*i_q + u_d/L_d, \
                    -i_q_dot-(L_d/L_q)*omega*i_d - (R_m/L_q)*i_q + u_q/L_q - (omega*K_m)/L_q)

    model = acados_dae()

    model.f_impl_expr = fimpl
    model.f_expl_expr = fexpl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = []
    model.p = p
    model.name = model_name

    return model 

# Hexagon Voltage Constraints
#====================================================================
def get_general_constraints_DC():

    # polytopic constraint on the input
    s3 = sqrt(3)
    cphi = cos(phi)
    sphi = sin(phi)
    
    h1 = s3*cphi + sphi
    h2 = sphi
    h3 = -s3*cphi + sphi
    h7 = -s3*sphi + cphi
    h8 = cphi
    h9 = s3*sphi + cphi

    # form D and C matrices
    # (acados C interface works with column major format)
    D = nmp.array([[h1, h7],[h2, h8],[h3, h9]])
    C = nmp.array([[0, 0],[0, 0],[0, 0]])

    g1 = 2.0/s3*u_max
    g2 = 1.0/s3*u_max

    lg = nmp.array([-g1, -g2, -g1])
    ug = nmp.array([g1, g2, g1])
 
    res = dict()
    res["D"] = D
    res["C"] = C
    res["lg"] = lg
    res["ug"] = ug

    return res

# Hexagon Terminal Constraints
#====================================================================
def get_general_terminal_constraints_DC():

    # polytopic constraint on the input
    s3 = sqrt(3)
    cphi = cos(phi)
    sphi = sin(phi)
    
    h1 = s3*cphi + sphi
    h2 = sphi
    h3 = -s3*cphi + sphi
    h7 = -s3*sphi + cphi
    h8 = cphi
    h9 = s3*sphi + cphi

    # form D and C matrices
    D = nmp.array([[h1, h7],[h2, h8],[h3, h9]])

    A = nmp.array([[-R_m/L_d, w_val*L_q/L_d],[-w_val*L_d/L_q, -R_m/L_q]])
    invB = nmp.array([[L_d,0],[0,L_q]])
    f = nmp.array([[0],[-K_m*w_val/L_q]])
    
    invBA = invB.dot(A)
    Ce = D.dot(invBA)

    g1 = 2.0/s3*u_max
    g2 = 1.0/s3*u_max

    invBf = invB.dot(f)
    lge = -nmp.array([[g1], [g2], [g1]]) - D.dot(invBf)
    uge = +nmp.array([[g1], [g2], [g1]]) - D.dot(invBf)
    
    res = dict()
    res["Ce"] = Ce
    res["lge"] = lge.flatten()
    res["uge"] = uge.flatten()

    return res

# # Torque Constraints
# #====================================================================
# def export_torqueline():

#     con_name = 'torqueline'

#     # set up states 
#     i_d = SX.sym('i_d')
#     i_q = SX.sym('i_q')
#     x = vertcat(i_d, i_q)

#     # set up controls
#     u_d = SX.sym('u_d')
#     u_q = SX.sym('u_q')
#     u = vertcat(u_d, u_q)

#     # set up parameters
#     omega  = SX.sym('omega')        # speed
#     dist_d = SX.sym('dist_d')       # d disturbance
#     dist_q = SX.sym('dist_q')       # q disturbance
#     tau_des = SX.sym('tau_des')     # Desired tau
#     p = vertcat(omega, dist_d, dist_q, tau_des)

#     # torque and voltage constraints
#     constraint = acados_constraint()

#     r = SX.sym('r', 1, 1)
#     constraint.con_h_expr = r[0]
#     constraint.con_r_expr = tau_des - 1.5*N_P*((L_d-L_q)*i_d*i_q + K_m*i_q)
#     constraint.x = x
#     constraint.u = u 
#     constraint.r = r
#     constraint.nr = 1
#     constraint.nh = 1
#     constraint.name = con_name

#     return constraint  

# Torque Constraints --> pd
#====================================================================
def export_torqueline_pd():

    con_name = 'torqueline'

    # set up states 
    i_d = SX.sym('i_d')
    i_q = SX.sym('i_q')
    x = vertcat(i_d, i_q)

    # set up controls
    u_d = SX.sym('u_d')
    u_q = SX.sym('u_q')
    u = vertcat(u_d, u_q)

    # set up parameters
    omega  = SX.sym('omega')        # speed
    dist_d = SX.sym('dist_d')       # d disturbance
    dist_q = SX.sym('dist_q')       # q disturbance
    tau_des = SX.sym('tau_des')     # Desired tau
    p = vertcat(omega, dist_d, dist_q, tau_des)

    # torque and voltage constraints
    constraint = acados_constraint()
    r = SX.sym('r', 3, 1)
    constraint.con_phi_expr = vertcat(r[0], r[1]**2 + r[2]**2)
    constraint.con_r_expr =  vertcat(tau_des - 1.5*N_P*((L_d-L_q)*i_d*i_q + K_m*i_q), u_d, u_q)
    constraint.x = x
    constraint.u = u 
    constraint.r = r
    constraint.p = p
    constraint.nr = 3
    constraint.nphi = 2
    constraint.name = con_name

    return constraint 

# # torque Constraints --> nl
# #====================================================================
# def export_torqueline_nl():

#     con_name = 'torqueline'

#     # set up states 
#     i_d = SX.sym('i_d')
#     i_q = SX.sym('i_q')
#     x = vertcat(i_d, i_q)

#     # set up controls
#     u_d = SX.sym('u_d')
#     u_q = SX.sym('u_q')
#     u = vertcat(u_d, u_q)

#     # set up parameters
#     omega  = SX.sym('omega')        # speed
#     dist_d = SX.sym('dist_d')       # d disturbance
#     dist_q = SX.sym('dist_q')       # q disturbance
#     tau_des = SX.sym('tau_des')     # Desired tau
#     p = vertcat(omega, dist_d, dist_q, tau_des)

#     # torque and voltage constraints
#     constraint = acados_constraint()
#     constraint.con_h_expr = tau_des - 1.5*N_P*((L_d-L_q)*i_d*i_q + K_m*i_q)
#     constraint.x = x
#     constraint.u = u 
#     constraint.nh = 1
#     constraint.name = con_name

#     return constraint  

# torque End-Constraints
#====================================================================
# def export_torquelineEnd():

#     con_name = 'torquelineEnd'

#     # set up states 
#     i_d = SX.sym('i_d')
#     i_q = SX.sym('i_q')
#     x = vertcat(i_d, i_q)

#     # set up controls
#     u_d = SX.sym('u_d')
#     u_q = SX.sym('u_q') 
#     u = vertcat(u_d, u_q)

#     # set up parameters
#     omega  = SX.sym('omega')        # speed
#     dist_d = SX.sym('dist_d')       # d disturbance
#     dist_q = SX.sym('dist_q')       # q disturbance
#     tau_des = SX.sym('tau_des')     # Desired tau
#     p = vertcat(omega, dist_d, dist_q, tau_des)
    
#     # torque and voltage constraints
#     constraint = acados_constraint()
#     r = SX.sym('r', 1, 1)
#     constraint.con_h_expr = r[0]
#     constraint.con_r_expr = tau_des - 1.5*N_P*((L_d-L_q)*i_d*i_q + K_m*i_q)
#     constraint.x = x
#     constraint.u = u 
#     constraint.r = r
#     constraint.nr = 1
#     constraint.nh = 1
#     constraint.name = con_name

#     return constraint

# torque End-Constraints --> pd
#====================================================================
def export_torquelineEnd_pd():

    con_name = 'torquelineEnd'

    # set up states 
    i_d = SX.sym('i_d')
    i_q = SX.sym('i_q')
    x = vertcat(i_d, i_q)

    # set up controls
    u_d = SX.sym('u_d')
    u_q = SX.sym('u_q') 
    u = vertcat(u_d, u_q)

    # set up parameters
    omega  = SX.sym('omega')        # speed
    dist_d = SX.sym('dist_d')       # d disturbance
    dist_q = SX.sym('dist_q')       # q disturbance
    tau_des = SX.sym('tau_des')     # Desired tau
    p = vertcat(omega, dist_d, dist_q, tau_des)
    
    alpha = R_m**2 + omega**2*L_d**2
    beta  = R_m**2 + omega**2*L_q**2
    gamma = 2*R_m*omega*(L_d-L_q)
    delta = 2*omega**2*L_d*K_m
    epsilon = 2*R_m*omega*K_m
    rho = omega**2*K_m**2

    # torque and voltage constraints
    constraint = acados_constraint()
    r = SX.sym('r', 3, 1)
    constraint.con_phi_expr = vertcat(r[0], alpha*r[1]**2 + beta*r[2]**2 + gamma*r[1]*r[2] + delta*r[1] + epsilon*r[2] + rho)
    constraint.con_r_expr = vertcat(tau_des - 1.5*N_P*((L_d-L_q)*i_d*i_q + K_m*i_q), i_d, i_q)
    constraint.x = x
    constraint.u = u 
    constraint.r = r
    constraint.p = p
    constraint.nr = 3
    constraint.nphi = 2
    constraint.name = con_name

    return constraint

# # torque End-Constraints --> nl
# #====================================================================
# def export_torquelineEnd_nl():

#     con_name = 'torquelineEnd'

#     # set up states 
#     i_d = SX.sym('i_d')
#     i_q = SX.sym('i_q')
#     x = vertcat(i_d, i_q)

#     # set up controls
#     u_d = SX.sym('u_d')
#     u_q = SX.sym('u_q') 
#     u = vertcat(u_d, u_q)

#     # set up parameters
#     omega  = SX.sym('omega')        # speed
#     dist_d = SX.sym('dist_d')       # d disturbance
#     dist_q = SX.sym('dist_q')       # q disturbance
#     tau_des = SX.sym('tau_des')     # Desired tau
#     p = vertcat(omega, dist_d, dist_q, tau_des)
    
#     # torque and voltage constraints
#     constraint = acados_constraint()
#     constraint.con_h_expr = tau_des - 1.5*N_P*((L_d-L_q)*i_d*i_q + K_m*i_q)
#     constraint.x = x
#     constraint.u = u
#     constraint.nh = 1
#     constraint.name = con_name

#     return constraint

# # |x - x_soll| Penalty
# #====================================================================
# def export_current_Penalty():

#     con_name = 'currentPenalty'

#     # set up states 
#     i_d = SX.sym('i_d')
#     i_q = SX.sym('i_q')
#     x = vertcat(i_d, i_q)

#     # set up controls
#     u_d = SX.sym('u_d')
#     u_q = SX.sym('u_q') 
#     u = vertcat(u_d, u_q)

#     # set up parameters
#     omega  = SX.sym('omega')        # speed
#     dist_d = SX.sym('dist_d')       # d disturbance
#     dist_q = SX.sym('dist_q')       # q disturbance
#     tau_des = SX.sym('tau_des')     # Desired tau
#     p = vertcat(omega, dist_d, dist_q, tau_des)
    
#     # torque and voltage constraints
#     constraint = acados_constraint()
#     constraint.con_h_expr = vertcat(i_d-i_d_ref, i_q-i_q_ref)
#     constraint.x = x
#     constraint.u = u
#     constraint.nh = 2
#     constraint.name = con_name

#     return constraint

# # |x - x_soll| Penalty -> pd
# #====================================================================
# def export_current_Penalty_pd():

#     con_name = 'currentPenalty'

#     # set up states 
#     i_d = SX.sym('i_d')
#     i_q = SX.sym('i_q')
#     x = vertcat(i_d, i_q)

#     # set up controls
#     u_d = SX.sym('u_d')
#     u_q = SX.sym('u_q')
#     u = vertcat(u_d, u_q)

#     # set up parameters
#     omega  = SX.sym('omega')        # speed
#     dist_d = SX.sym('dist_d')       # d disturbance
#     dist_q = SX.sym('dist_q')       # q disturbance
#     tau_des = SX.sym('tau_des')     # Desired tau
#     p = vertcat(omega, dist_d, dist_q, tau_des)

#     # torque and voltage constraints
#     constraint = acados_constraint()
#     r = SX.sym('r', 4, 1)
#     constraint.con_h_expr = vertcat(r[0], r[1], r[2]**2 + r[3]**2)
#     constraint.con_r_expr =  vertcat(i_d-i_d_ref, i_q-i_q_ref, u_d, u_q)
#     constraint.x = x
#     constraint.u = u 
#     constraint.r = r
#     constraint.nr = 4
#     constraint.nh = 3
#     constraint.name = con_name

#     return constraint 

# # |x - x_soll| End Penalty -> pd
# #====================================================================
# def export_current_PenaltyEnd_pd():

#     con_name = 'currentPenaltyEnd'

#     # set up states 
#     i_d = SX.sym('i_d')
#     i_q = SX.sym('i_q')
#     x = vertcat(i_d, i_q)

#     # set up controls
#     u_d = SX.sym('u_d')
#     u_q = SX.sym('u_q') 
#     u = vertcat(u_d, u_q)

#     # set up parameters
#     omega  = SX.sym('omega')        # speed
#     dist_d = SX.sym('dist_d')       # d disturbance
#     dist_q = SX.sym('dist_q')       # q disturbance
#     tau_des = SX.sym('tau_des')     # Desired tau
#     p = vertcat(omega, dist_d, dist_q, tau_des)
    
#     # torque and voltage constraints
#     constraint = acados_constraint()
#     r = SX.sym('r', 2, 1)
#     constraint.con_h_expr = vertcat(r[0], r[1])
#     constraint.con_r_expr = vertcat(i_d-i_d_ref, i_q-i_q_ref)
#     constraint.x = x
#     constraint.u = u 
#     constraint.r = r
#     constraint.nr = 2
#     constraint.nh = 2
#     constraint.name = con_name

#     return constraint

# create render arguments
ra = acados_ocp_nlp()

# export model
model = export_ode_model()

# set model_name
ra.model = model

# export constraint description
# if FORMULATION == 1:
#     constraint_nltorqueline = export_torqueline()
#     constraint_nltorquelineEnd = export_torquelineEnd()
# if FORMULATION == 2:
#     constraint_nltorqueline = export_torqueline_pd()
#     constraint_nltorquelineEnd = export_torquelineEnd()
# if FORMULATION == 3:
#     constraint_nltorqueline = export_torqueline_nl()
#     constraint_nltorquelineEnd = export_torquelineEnd_nl()
if FORMULATION == 4:
    constraint_nltorqueline = export_torqueline_pd()
    constraint_nltorquelineEnd = export_torquelineEnd_pd()
# if FORMULATION == 5:
#     constraint_nltorqueline = export_current_Penalty()
#     constraint_nltorquelineEnd = export_current_Penalty()
# if FORMULATION == 6:
#     constraint_nltorqueline = export_torqueline_pd()
#     constraint_nltorquelineEnd = export_torquelineEnd_nl()
# if FORMULATION == 7:
#     constraint_nltorqueline = export_current_Penalty_pd()
#     constraint_nltorquelineEnd = export_current_PenaltyEnd_pd()

# set model dims
nx = model.x.size()[0]
nu = model.u.size()[0]
np = model.p.size()[0]
ny = nu + nx
ny_e = nx
Tf = N*Ts

# set ocp_nlp_dimensions
nlp_dims = ra.dims
nlp_dims.nx = nx
nlp_dims.nu = nu
nlp_dims.np = np
nlp_dims.ny = ny
nlp_dims.ny_e = ny_e

if FORMULATION == 0:
    nlp_dims.nbu = 0
    nlp_dims.ng = 3
    nlp_dims.ng_e = 0

# if FORMULATION == 1:        
#     nlp_dims.nbu = 0
#     nlp_dims.ng = 3
#     nlp_dims.ng_e = 3
#     nlp_dims.nh_e = 1                   # 1st torque constraint 
#     nlp_dims.npd_e = 1                  # 1st torque constraint 
#     nlp_dims.ns = 1
#     nlp_dims.ns_e = 1 
#     nlp_dims.nsh = 1
#     nlp_dims.nsh_e = 1 

# if FORMULATION == 2:        
#     nlp_dims.nbu = 0
#     nlp_dims.ng = 3
#     nlp_dims.ng_e = 3
#     nlp_dims.nh  = 2            # 1st torque constraint | 2nd voltage constraint
#     nlp_dims.npd  = 3           # 1st torque constraint | 2nd voltage constraint
#     nlp_dims.nh_e = 1           # 1st torque constraint | 2nd terminal set
#     nlp_dims.npd_e = 0          # current constraint
#     nlp_dims.ns = 2
#     nlp_dims.ns_e = 2 
#     nlp_dims.nsh = 2
#     nlp_dims.nsh_e = 2 

if FORMULATION == 4:        
    nlp_dims.nbu = 0
    nlp_dims.ng = 3
    nlp_dims.ng_e = 3
    nlp_dims.nphi  = 2            # 1st torque constraint | 2nd voltage constraint
    nlp_dims.nr = 3           # 1st torque constraint | 2nd voltage constraint
    nlp_dims.nphi_e = 2           # 1st torque constraint | 2nd terminal set
    nlp_dims.nr_e = 3          # 1st torque constraint | 2nd terminal set
    nlp_dims.ns = 1
    nlp_dims.ns_e = 1 
    nlp_dims.nsphi = 1
    nlp_dims.nsphi_e = 1 

# if FORMULATION == 6:        
#     nlp_dims.nbu = 0
#     nlp_dims.ng = 3
#     nlp_dims.ng_e = 3
#     nlp_dims.nh  = 2            # 1st torque constraint | voltage constraint
#     nlp_dims.npd  = 3           # 1st torque constraint | voltage constraint
#     nlp_dims.nh_e = 1           # 1st torque constraint 
#     nlp_dims.npd_e = 0          # 1st torque constraint
#     nlp_dims.ns = 1
#     nlp_dims.ns_e = 1 
#     nlp_dims.nsh = 1
#     nlp_dims.nsh_e = 1 

# if FORMULATION == 7:        
#     nlp_dims.nbu = 0
#     nlp_dims.ng = 3
#     nlp_dims.ng_e = 0
#     nlp_dims.nh  = 3            # current constraint | voltage constraint
#     nlp_dims.npd  = 4           # current constraint | voltage constraint
#     nlp_dims.nh_e = 2           # current constraint
#     nlp_dims.npd_e = 2          # current constraint
#     nlp_dims.ns = 2
#     nlp_dims.ns_e = 2 
#     nlp_dims.nsh = 2
#     nlp_dims.nsh_e = 2 

nlp_dims.nbx = 0
nlp_dims.nbx_e = 0
nlp_dims.N = N

# set weighting matrices
nlp_cost = ra.cost

Q = nmp.eye(nx)
Q[0,0] = wd*Weight_TUNING*Tf/N
Q[1,1] = wq*Weight_TUNING*Tf/N

Q_e = nmp.eye(nx)
Q_e[0,0] = wd*Weight_E_TUNING*Tf/N
Q_e[1,1] = wq*Weight_E_TUNING*Tf/N

R = nmp.eye(nu)
R[0,0] = INPUT_REG*Tf/N
R[1,1] = INPUT_REG*Tf/N

nlp_cost.W = scipy.linalg.block_diag(Q, R)   # weight matrix
nlp_cost.W_e = Q_e                           # weight matrix for Mayer term

Vu = nmp.zeros((ny, nu))
Vu[2,0] = 1.0
Vu[3,1] = 1.0
nlp_cost.Vu = Vu

Vx = nmp.zeros((ny, nx))
Vx[0,0] = 1.0
Vx[1,1] = 1.0
nlp_cost.Vx = Vx

Vx_e = nmp.zeros((ny_e, nx))
Vx_e[0,0] = 1.0
Vx_e[1,1] = 1.0
nlp_cost.Vx_e = Vx_e

if FORMULATION == 0 or FORMULATION == 5 or FORMULATION == 6 or FORMULATION == 7:
    nlp_cost.yref  = nmp.zeros((ny, ))
    nlp_cost.yref[0] = i_d_ref
    nlp_cost.yref[1] = i_q_ref
    nlp_cost.yref[2] = 0
    nlp_cost.yref[3] = 0
    nlp_cost.yref_e = nmp.zeros((ny_e, ))
    nlp_cost.yref_e[0] = i_d_ref
    nlp_cost.yref_e[1] = i_q_ref
else: 
    nlp_cost.yref  = nmp.zeros((ny, ))
    nlp_cost.yref[0] = 0
    nlp_cost.yref[1] = 0
    nlp_cost.yref[2] = 0
    nlp_cost.yref[3] = 0
    nlp_cost.yref_e = nmp.zeros((ny_e, ))
    nlp_cost.yref_e[0] = 0
    nlp_cost.yref_e[1] = 0

# get D and C
res = get_general_constraints_DC()
D = res["D"]
C = res["C"]
lg = res["lg"]
ug = res["ug"]

res = get_general_terminal_constraints_DC()
Ce = res["Ce"]
lge = res["lge"]
uge = res["uge"]

# setting bounds
nlp_con = ra.constraints

if FORMULATION == 0: 
    # setting general constraints --> lg <= D*u + C*u <= ug
    nlp_con.D   = D
    nlp_con.C   = C
    nlp_con.lg  = lg
    nlp_con.ug  = ug

    nlp_con.x0 = x0Start

# if FORMULATION == 1 or FORMULATION == 3: 
#     # setting general constraints --> lg <= D*u + C*u <= ug
#     nlp_con.D   = D
#     nlp_con.C   = C
#     nlp_con.lg  = lg
#     nlp_con.ug  = ug
#     nlp_con.C_e  = Ce
#     nlp_con.lg_e = lge
#     nlp_con.ug_e = uge

#     # lower gradient/hessian wrt lower slack
#     nlp_cost.zl = SLACK_TUNING*nmp.array([1])
#     nlp_cost.Zl = SLACK_TUNINGHessian*nmp.ones((1,1))                # hessian
#     nlp_cost.zu = SLACK_TUNING*nmp.array([1])
#     nlp_cost.Zu = SLACK_TUNINGHessian*nmp.ones((1,1))                # hessian
#     #_e
#     nlp_cost.zl_e = SLACK_E_TUNING*nmp.array([1])           
#     nlp_cost.Zl_e = SLACK_TUNINGHessian*nmp.ones((1,1))              # hessian   
#     nlp_cost.zu_e = SLACK_E_TUNING*nmp.array([1])     
#     nlp_cost.Zu_e = SLACK_TUNINGHessian*nmp.ones((1,1))              # hessian   

#     nlp_con.lh = nmp.array([0])                         # 1st torque constraint 
#     nlp_con.uh = nmp.array([0])
#     nlp_con.lsh = nmp.array([0])                        # soft lower bounds --> greater than 0
#     nlp_con.ush = nmp.array([0])                        # soft upper bounds --> greater than 0
#     #_e
#     nlp_con.lh_e = nmp.array([0])                       # 1st torque constraint 
#     nlp_con.uh_e = nmp.array([0])
#     nlp_con.lsh_e = nmp.array([0])      
#     nlp_con.ush_e = nmp.array([0])      

#     nlp_con.idxsh = nmp.array([0])
#     nlp_con.idxsh_e = nmp.array([0])  

#     nlp_con.x0 = x0Start

# if FORMULATION == 2: 
#     # setting general constraints --> lg <= D*u + C*u <= ug
#     nlp_con.D   = D
#     nlp_con.C   = C
#     nlp_con.lg  = lg
#     nlp_con.ug  = ug
#     nlp_con.C_e  = Ce
#     nlp_con.lg_e = lge
#     nlp_con.ug_e = uge

#     # lower gradient/hessian wrt lower slack
#     nlp_cost.zl = SLACK_TUNING*nmp.array([1])
#     nlp_cost.Zl = SLACK_TUNINGHessian*nmp.ones((1,1))                # hessian
#     nlp_cost.zu = SLACK_TUNING*nmp.array([1])
#     nlp_cost.Zu = SLACK_TUNINGHessian*nmp.ones((1,1))                # hessian
#     #_e
#     nlp_cost.zl_e = SLACK_E_TUNING*nmp.array([1])           
#     nlp_cost.Zl_e = SLACK_TUNINGHessian*nmp.ones((1,1))              # hessian   
#     nlp_cost.zu_e = SLACK_E_TUNING*nmp.array([1])     
#     nlp_cost.Zu_e = SLACK_TUNINGHessian*nmp.ones((1,1))              # hessian   

#     nlp_con.lh = nmp.array([0, -1e9])                   # 1st torque constraint | 2nd voltage constraint 
#     nlp_con.uh = nmp.array([0, u_max**2/3])
#     nlp_con.lsh = nmp.array([0])                        # soft lower bounds --> greater than 0
#     nlp_con.ush = nmp.array([0])                        # soft upper bounds --> greater than 0
#     #_e
#     nlp_con.lh_e = nmp.array([0])                       # 1st torque constraint | 2nd terminal set 
#     nlp_con.uh_e = nmp.array([0])
#     nlp_con.lsh_e = nmp.array([0])      
#     nlp_con.ush_e = nmp.array([0])      

#     nlp_con.idxsh = nmp.array([0])
#     nlp_con.idxsh_e = nmp.array([0])  

#     nlp_con.x0 = x0Start

if FORMULATION == 4: 
    # setting general constraints --> lg <= D*u + C*u <= ug
    nlp_con.D   = D
    nlp_con.C   = C
    nlp_con.lg  = lg
    nlp_con.ug  = ug
    nlp_con.C_e  = Ce
    nlp_con.lg_e = lge
    nlp_con.ug_e = uge

    # lower gradient/hessian wrt lower slack
    nlp_cost.zl = SLACK_TUNING*nmp.array([1])
    nlp_cost.Zl = SLACK_TUNINGHessian*nmp.ones((1,))                # hessian
    nlp_cost.zu = SLACK_TUNING*nmp.array([1])
    nlp_cost.Zu = SLACK_TUNINGHessian*nmp.ones((1,))                # hessian
    #_e
    nlp_cost.zl_e = SLACK_E_TUNING*nmp.array([1])           
    nlp_cost.Zl_e = SLACK_TUNINGHessian*nmp.ones((1,))              # hessian   
    nlp_cost.zu_e = SLACK_E_TUNING*nmp.array([1])     
    nlp_cost.Zu_e = SLACK_TUNINGHessian*nmp.ones((1,))              # hessian   

    nlp_con.lphi = nmp.array([0, -1e9])                       # 1st torque constraint | 2nd voltage constraint 
    nlp_con.uphi = nmp.array([0, u_max**2/3])
    nlp_con.lsphi = nmp.array([0])                            # soft lower bounds --> greater than 0
    nlp_con.usphi = nmp.array([0])                            # soft upper bounds --> greater than 0
    #_e
    nlp_con.lphi_e = nmp.array([0, -1e9])                     # 1st torque constraint | 2nd terminal set 
    nlp_con.uphi_e = nmp.array([0, u_max**2/3])
    nlp_con.lsphi_e = nmp.array([0])      
    nlp_con.usphi_e = nmp.array([0])      

    nlp_con.idxsphi = nmp.array([0])
    nlp_con.idxsphi_e = nmp.array([0])  

    nlp_con.x0 = x0Start

# if FORMULATION == 5:
#     # setting general constraints --> lg <= D*u + C*u <= ug
#     nlp_con.D   = D
#     nlp_con.C   = C
#     nlp_con.lg  = lg
#     nlp_con.ug  = ug
#     nlp_con.C_e  = Ce
#     nlp_con.lg_e = lge
#     nlp_con.ug_e = uge

#     # lower gradient/hessian wrt lower slack
#     nlp_cost.zl = SLACK_TUNING*nmp.array([1,1])
#     nlp_cost.Zl = SLACK_TUNINGHessian*nmp.ones((2,2))                # hessian
#     nlp_cost.zu = SLACK_TUNING*nmp.array([1,1])
#     nlp_cost.Zu = SLACK_TUNINGHessian*nmp.ones((2,2))                # hessian
#     #_e
#     nlp_cost.zl_e = SLACK_E_TUNING*nmp.array([1,1])           
#     nlp_cost.Zl_e = SLACK_TUNINGHessian*nmp.ones((2,2))              # hessian   
#     nlp_cost.zu_e = SLACK_E_TUNING*nmp.array([1,1])     
#     nlp_cost.Zu_e = SLACK_TUNINGHessian*nmp.ones((2,2))              # hessian   

#     nlp_con.lh = nmp.array([0,0])                     # 1st current constraint  
#     nlp_con.uh = nmp.array([0,0])
#     nlp_con.lsh = nmp.array([0,0])                    # soft lower bounds --> greater than 0
#     nlp_con.ush = nmp.array([0,0])                    # soft upper bounds --> greater than 0
#     #_e
#     nlp_con.lh_e = nmp.array([0,0])                   # 1st current constraint  
#     nlp_con.uh_e = nmp.array([0,0])
#     nlp_con.lsh_e = nmp.array([0,0])      
#     nlp_con.ush_e = nmp.array([0,0])      

#     nlp_con.idxsh = nmp.array([0,1])
#     nlp_con.idxsh_e = nmp.array([0,1])  

#     nlp_con.x0 = x0Start

# if FORMULATION == 6: 
#     # setting general constraints --> lg <= D*u + C*u <= ug
#     nlp_con.D   = D
#     nlp_con.C   = C
#     nlp_con.lg  = lg
#     nlp_con.ug  = ug
#     nlp_con.C_e  = Ce
#     nlp_con.lg_e = lge
#     nlp_con.ug_e = uge

#     # lower gradient/hessian wrt lower slack
#     nlp_cost.zl = SLACK_TUNING*nmp.array([1])
#     nlp_cost.Zl = SLACK_TUNINGHessian*nmp.ones((1,1))                # hessian
#     nlp_cost.zu = SLACK_TUNING*nmp.array([1])
#     nlp_cost.Zu = SLACK_TUNINGHessian*nmp.ones((1,1))                # hessian
#     #_e
#     nlp_cost.zl_e = SLACK_E_TUNING*nmp.array([1])           
#     nlp_cost.Zl_e = SLACK_TUNINGHessian*nmp.ones((1,1))              # hessian   
#     nlp_cost.zu_e = SLACK_E_TUNING*nmp.array([1])     
#     nlp_cost.Zu_e = SLACK_TUNINGHessian*nmp.ones((1,1))              # hessian   

#     nlp_con.lh = nmp.array([0, -1e9])                   # 1st torque constraint | 2nd voltage constraint 
#     nlp_con.uh = nmp.array([0, u_max**2/3])
#     nlp_con.lsh = nmp.array([0])                        # soft lower bounds --> greater than 0
#     nlp_con.ush = nmp.array([0])                        # soft upper bounds --> greater than 0
#     #_e
#     nlp_con.lh_e = nmp.array([0])                       # 1st torque constraint
#     nlp_con.uh_e = nmp.array([0])
#     nlp_con.lsh_e = nmp.array([0])      
#     nlp_con.ush_e = nmp.array([0])      

#     nlp_con.idxsh = nmp.array([0])
#     nlp_con.idxsh_e = nmp.array([0])  

#     nlp_con.x0 = x0Start

# if FORMULATION == 7:
#     # setting general constraints --> lg <= D*u + C*u <= ug
#     nlp_con.D   = D
#     nlp_con.C   = C
#     nlp_con.lg  = lg
#     nlp_con.ug  = ug
#     # nlp_con.C_e  = Ce
#     # nlp_con.lg_e = lge
#     # nlp_con.ug_e = uge

#     # lower gradient/hessian wrt lower slack
#     nlp_cost.zl = SLACK_TUNING*nmp.array([1,1])
#     nlp_cost.Zl = SLACK_TUNINGHessian*nmp.ones((2,2))                # hessian
#     nlp_cost.zu = SLACK_TUNING*nmp.array([1,1])
#     nlp_cost.Zu = SLACK_TUNINGHessian*nmp.ones((2,2))                # hessian
#     #_e
#     nlp_cost.zl_e = SLACK_E_TUNING*nmp.array([1,1])           
#     nlp_cost.Zl_e = SLACK_TUNINGHessian*nmp.ones((2,2))              # hessian   
#     nlp_cost.zu_e = SLACK_E_TUNING*nmp.array([1,1])     
#     nlp_cost.Zu_e = SLACK_TUNINGHessian*nmp.ones((2,2))              # hessian   

#     nlp_con.lh = nmp.array([0,0,-1e9])                      # 1st current constraint | 2nd voltage constraint  
#     nlp_con.uh = nmp.array([0,0,u_max**2/3])
#     nlp_con.lsh = nmp.array([0,0])                          # soft lower bounds --> greater than 0
#     nlp_con.ush = nmp.array([0,0])                          # soft upper bounds --> greater than 0
#     #_e
#     nlp_con.lh_e = nmp.array([0,0])                         # 1st current constraint  
#     nlp_con.uh_e = nmp.array([0,0])
#     nlp_con.lsh_e = nmp.array([0,0])      
#     nlp_con.ush_e = nmp.array([0,0])      

#     nlp_con.idxsh = nmp.array([0,1])
#     nlp_con.idxsh_e = nmp.array([0,1])  

#     nlp_con.x0 = x0Start

# setting parameters
nlp_con.p = nmp.array([w_val, 0.0, 0.0, tau_wal])

# set QP solver
# ra.solver_config.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
ra.solver_config.qp_solver = 'FULL_CONDENSING_HPIPM'
# ra.solver_config.qp_solver = 'FULL_CONDENSING_QPOASES'
ra.solver_config.hessian_approx = 'GAUSS_NEWTON'
ra.solver_config.integrator_type = 'IRK'
# ra.solver_config.integrator_type = 'ERK'
ra.solver_config.sim_method_num_stages = 1  # 1: RK1, 2: RK2, 4: RK4    

# set prediction horizon
ra.solver_config.tf = Tf
ra.solver_config.nlp_solver_type = 'SQP_RTI'
# ra.solver_config.nlp_solver_type = 'SQP'

# set header path
# ra.acados_include_path = '/mnt/c/Austauch/projects/SCQPAndrea/acados/include'
# ra.acados_lib_path = '/mnt/c/Austauch/projects/SCQPAndrea/acados/lib'
ra.acados_include_path = '~/acados/include'
ra.acados_lib_path = '~/acados/lib'

file_name = 'acados_ocp.json'

# import pdb; pdb.set_trace()

if CODE_GEN == 1:
    if FORMULATION == 0:
        acados_solver = generate_solver(ra, json_file = file_name)
    # if FORMULATION == 1:
    #     ra.con_h = constraint_nltorqueline
    #     ra.con_h_e = constraint_nltorquelineEnd
    #     acados_solver = generate_solver(ra, json_file = file_name)
    # if FORMULATION == 2:
    #     ra.con_h = constraint_nltorqueline
    #     ra.con_h_e = constraint_nltorquelineEnd
    #     acados_solver = generate_solver(ra, json_file = file_name)
    # if FORMULATION == 3:
    #     ra.con_h = constraint_nltorqueline
    #     ra.con_h_e = constraint_nltorquelineEnd
    #     acados_solver = generate_solver(ra, json_file = file_name)
    if FORMULATION == 4:
        ra.con_phi = constraint_nltorqueline
        ra.con_phi_e = constraint_nltorquelineEnd
        nlp_con.constr_type = 'BGP'
        nlp_con.constr_type_e = 'BGP'
        acados_solver = generate_solver(ra, json_file = file_name)
    # if FORMULATION == 5:
    #     ra.con_h = constraint_nltorqueline
    #     ra.con_h_e = constraint_nltorquelineEnd
    #     acados_solver = generate_solver(ra, json_file = file_name)
    # if FORMULATION == 6:
    #     ra.con_h = constraint_nltorqueline
    #     ra.con_h_e = constraint_nltorquelineEnd
    #     acados_solver = generate_solver(ra, json_file = file_name)
    # if FORMULATION == 7:
    #     ra.con_h = constraint_nltorqueline
    #     ra.con_h_e = constraint_nltorquelineEnd
    #     acados_solver = generate_solver(ra, json_file = file_name)

if COMPILE == 1:
    # make 
    os.chdir('c_generated_code')
    os.system('make')
    os.system('make shared_lib')
    os.chdir('..')

# closed loop simulation
Nsim = 20

simX = nmp.ndarray((Nsim, nx))
simU = nmp.ndarray((Nsim, nu))
simXR = nmp.ndarray((Nsim+1, nx))
simXRN = nmp.ndarray((Nsim, nx))

print("============================================================================================")
print("Mode = ", FORMULATION)
print("speed_el = ", w_val)
print("Sample Time = ", Ts)
print("============================================================================================")
print("Initial Condition")
print("id0 = ", x0Start[0])
print("iq0 = ", x0Start[1])

# get initial condition for real simulation
simXR[0,0] = x0Start[0]
simXR[0,1] = x0Start[1]

xvec = nmp.matrix([[x0Start[0]],[x0Start[1]]])

# compute warm-start
for i in range(WARMSTART_ITERS):
    status = acados_solver.solve()

for i in range(Nsim):
    print("\n")
    print("SimulationStep = ", i)
    print("=================")
    status = acados_solver.solve()

    if status != 0:
        print('acados returned status {}. Exiting.'.format(status))
        import pdb; pdb.set_trace()

    # get solution
    x0 = acados_solver.get(0, "x")
    xN = acados_solver.get(N, "x")
    u0 = acados_solver.get(0, "u")
    
    for j in range(nx):
        simX[i,j] = x0[j]
        simXRN[i,j] = xN[j]

    for j in range(nu):
        simU[i,j] = u0[j]

    uvec = nmp.matrix([[u0[0]],[u0[1]]])

    # real Simulation 
    A = nmp.matrix([[-R_m/L_d, w_val*L_q/L_d],[-w_val*L_d/L_q, -R_m/L_q]],dtype=float)
    B = nmp.matrix([[1.0/L_d,0.0],[0.0,1.0/L_q]],dtype=float)
    f = nmp.matrix([[0.0],[-K_m*w_val/L_q]],dtype=float)

    # Euler
    #========================
    # xvec = xvec + Ts*A*xvec + Ts*B*uvec + Ts*f

    # Heun
    #========================
    # xvec = xvec + Ts*A*xvec + Ts*B*uvec + Ts*f + 0.5*Ts*Ts*A*A*xvec + 0.5*Ts*Ts*A*B*uvec + 0.5*Ts*Ts*A*f

    # Z-Transformation
    #========================
    Ad = scipy.linalg.expm(A*Ts_sim)
    invA = nmp.linalg.inv(A)
    Bd = invA*(Ad-nmp.eye(2))*B
    fd = invA*(Ad-nmp.eye(2))*f
    xvec = Ad*xvec + Bd*uvec + fd   

    print("States= ", xvec)
    print("Controls= ", uvec)
    print("\n")

    # update initial condition xk+1
    acados_solver.set(0, "lbx",  xvec)
    acados_solver.set(0, "ubx",  xvec)

    # update initial condition
    # ============================================
    # x0 = acados_solver.get(1, "x")
    # acados_solver.set(0, "lbx",  x0)
    # acados_solver.set(0, "ubx",  x0)

    simXR[i+1,0] = xvec[0]
    simXR[i+1,1] = xvec[1]

# plot results
t = nmp.linspace(0.0, Ts*Nsim, Nsim)
plt.subplot(4, 1, 1)
plt.step(t, simU[:,0], 'r')
plt.step(t, simU[:,0], 'ro')
plt.plot([0, Ts*Nsim], [nlp_cost.yref[2], nlp_cost.yref[2]], '--')
plt.title('closed-loop simulation')
plt.ylabel('u_d')
plt.xlabel('t')
plt.grid(True)
plt.subplot(4, 1, 2)
plt.step(t, simU[:,1], 'r')
plt.step(t, simU[:,1], 'ro')
plt.plot([0, Ts*Nsim], [nlp_cost.yref[3], nlp_cost.yref[3]], '--')
plt.ylabel('u_q')
plt.xlabel('t')
plt.grid(True)
plt.subplot(4, 1, 3)
plt.plot(t, simX[:,0])
plt.ylabel('x_d')
plt.xlabel('t')
plt.grid(True)
plt.subplot(4, 1, 4)
plt.plot(t, simX[:,1])
plt.ylabel('x_q')
plt.xlabel('t')
plt.grid(True)

# plot hexagon
r = 2/3*u_max
x1 = r
y1 = 0
x2 = r*cos(pi/3)
y2 = r*sin(pi/3)
q1 = -(y2 - y1/x1*x2)/(1-x2/x1)
m1 = -(y1 + q1)/x1

# box constraints
m2 = 0
q2 = r*sin(pi/3)
# -q2 <= uq  <= q2

plt.figure()
plt.plot(simU[:,0], simU[:,1], 'o')
plt.xlabel('ud')
plt.ylabel('uq')
ud = nmp.linspace(-1.5*u_max, 1.5*u_max, 100)
plt.plot(ud, -m1*ud -q1)
plt.plot(ud, -m1*ud +q1)
plt.plot(ud, +m1*ud -q1)
plt.plot(ud, +m1*ud +q1)
plt.plot(ud, -q2*nmp.ones((100, 1)))
plt.plot(ud, q2*nmp.ones((100, 1)))
plt.grid(True)
ax = plt.gca()
ax.set_xlim([-u_max, u_max])
ax.set_ylim([-u_max, u_max])
circle = plt.Circle((0, 0), u_max/nmp.sqrt(3), color='red', fill=False)
ax.add_artist(circle)

delta = 100
x = nmp.linspace(-i_max, i_max/3, delta)
y = nmp.linspace(-i_max, i_max, delta)
XV, YV = nmp.meshgrid(x,y)

alpha = R_m**2 + w_val**2*L_d**2
beta  = R_m**2 + w_val**2*L_q**2
gamma = 2*R_m*w_val*(L_d-L_q)
delta = 2*w_val**2*L_d*K_m
epsilon = 2*R_m*w_val*K_m
rho = w_val**2*K_m**2 - (u_max)**2/3

FeaSetV = alpha*XV**2 + beta*YV**2 + gamma*XV*YV + delta*XV + epsilon*YV + rho
TauV = 1.5*N_P*((L_d-L_q)*XV*YV + K_m*YV)

# trajectory in the dq-plane
plt.figure()
plt.plot(simXR[:,0], simXR[:,1],'bo')
plt.plot(simXR[:,0], simXR[:,1],'b')
plt.plot(simX[:,0], simX[:,1],'r')
plt.plot(simX[:,0], simX[:,1],'r*')
plt.contour(XV, YV, FeaSetV,[0],colors='r')
plt.contour(XV, YV, TauV,[tau_wal],colors='g')
plt.xlabel('x_d')
plt.ylabel('x_q')
plt.grid(True)

S1 = -0.24543672*XV + 0.50146524*YV
S2 = -0.214*XV      - 0.01815*YV
S3 = -0.18256328*XV -0.53776524*YV

plt.figure()
plt.plot(simXRN[:,0], simXRN[:,1],'bo')
plt.plot(simXRN[:,0], simXRN[:,1],'b')
plt.contour(XV, YV, FeaSetV,[0],colors='r')
plt.contour(XV, YV, S1,[-27.82562584,83.02562584],colors='k')
plt.contour(XV, YV, S2,[-0.11281292,55.31281292], colors='k')
plt.contour(XV, YV, S3,[-27.82562584,83.02562584],colors='k')
plt.contour(XV, YV, TauV,[tau_wal],colors='g')
plt.xlabel('x_d')
plt.ylabel('x_q')
plt.grid(True)

plt.show()

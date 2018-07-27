
from casadi import *

x, u, z = SX.sym('x', 2), SX.sym('u', 2), SX.sym('z', 2)
u_r, slack = SX.sym('u', 2), SX.sym('slack')
d1 = SX.sym('speed')

states = vertcat(u, x)
u1 = states[0]
u2 = states[1]
xD1_scaled = states[2]
xD2_scaled = states[3]

controls = vertcat(u_r, slack)
u1_r = controls[0]
u2_r = controls[1]

xA1 = z[0]
xA2 = z[1]

nx = states.numel()
nu = controls.numel()

Ts = 0.1

# Konstanten
m_dot_gas = 0               # M DOT GAS
T_ut = 700                  # Temperature Upstream Turbine
T_ima_k = 300               # Temperature IMA
m_dot_diesel = 0  
p_dt = 1                    # 1e5
# gamma=0.3904
p0 = 1                      # 1e5
# Rs_air = 287.1200         # gas constant air
# Rs_methane = 518.3600     # gas constant methane
# Vd = 4.9200e-04           # displacement volume 1 cylinder
# Vc = 3.1742e-05           # crevice volume
# gamma0 = 0.9070           # parameters vol efficiency speed dependency
# gamma1 = -5.6570e-05      # parameters vol efficiency speed dependency
kappa = 1.3                 # 1.4       # heat capacity ratio
c_p = 1050
eta_t = 0.5
eta_c = 0.545
T_amb = 20 + 273.15
J = 2e-7                    # 5e-7 # 1e-7 # 2e10
EPS = 1e-4
M_l = 28.97                 # g/mol
R_l = 8.314/M_l * 1e3       # J/(kg K)
R_air = 287.12
V_ident = 0.0001            # 0.001
# M_air = 28.9697
# M_exh = M_air
# M_O2 = 32.00
# XO2air = 0.21
# R_gas_stoich=17.16
# R_diesel_stoich=14.5
    
# Eingaenge skalieren, Zwischengroessen berechnen und ggf. limiteren
# Turbolader-Drehzahl [rpm] - Berechnung und Begrenzung
n_tc1 = xD1_scaled*1e5      # Turboladerdrehzahl in 1/min
n_tc_sqr = n_tc1**2
n_tc_sqr_scale = 1/2*sqrt(n_tc_sqr**2 + EPS**2) + 1/2*n_tc_sqr
n_tc = 1/2*sqrt((n_tc1-0.1*1e5)**2 + EPS)  + 1/2*(n_tc1-0.1*1e5)+0.1*1e5
# Einlassdruck-Berechnung
p_im = n_tc_sqr_scale/1e10*0.3811+0.944
# Motor-Drehzahl Umrechnung w_e [rad/s] -> speed_rpm [rpm]
w_e = d1
speed_rpm= w_e /(2*pi)*60;  # wird nicht benutzt...
# Druckverhaeltnisse -  Berechnung und Begrenzung
p_em = xD2_scaled
pi_m = p_em/p_im            # Druckverhaeltnis ueber Motor
pi_m_s = 1/2*sqrt(pi_m**2 + EPS) + 1/2*pi_m
pi_t = p_em/p_dt            # Druckverhaeltnis ueber Turbine
pi_t_s = 1/2*sqrt((pi_t-1.01)**2 + EPS) + 1/2*(pi_t-1.01)+1.01
    
# m dot beta und Exhaust
rho_im = p_im*1e5/(R_air*T_ima_k)
# deviation from nominal engine speed (2000rpm)
dz = w_e-2000*pi/30
mdot_beta = (0.005145) + (0.011)*rho_im + (0.01099)*pi_m_s + (1.1974e-04)*dz

m_dot_exh = m_dot_gas+m_dot_diesel+mdot_beta
m_dot_egr_hp = 0
    
# Massenfluss Turbine
# u_VTG = 1/2*sqrt((u1).**2 + EPS) + 1/2.*(u1);
u_VTG_lim = 1/2*sqrt(u1**2 + EPS) + 1/2*u1
u_VTG_shift = 1/2*sqrt((u1-1)**2 + EPS) + 1/2*(u1-1)
u_VTG = u_VTG_lim - u_VTG_shift
    
p00 =     -0.1561
p10 =      0.2414
p01 =      0.2819
p20 =     -0.2455
p11 =     -0.2829
p02 =     -0.1244
p30 =      0.1076
p21 =       0.199
p12 =     0.05707
p03 =     0.02601
p31 =    -0.08487
p22 =   -0.000165
p13 =   -0.007593
p04 =   -0.001824

mdot_turbCorr = p00 + p10*u_VTG + p01*pi_t_s + p20*u_VTG**2 + p11*u_VTG*pi_t_s + p02*pi_t_s**2 \
                    + p30*u_VTG**3 + p21*u_VTG**2*pi_t_s + p12*u_VTG*pi_t_s**2 + p03*pi_t_s**3 \
                    + p31*u_VTG**3*pi_t_s + p22*u_VTG**2*pi_t_s**2 + p13*u_VTG*pi_t_s**3 + p04*pi_t_s**4

mdot_turb=mdot_turbCorr*(p_em/p0)*sqrt(T_amb/T_ut)
    

# TC Leistungen     
P_t = mdot_turb*c_p*T_ut*eta_t*(1-pi_t_s**((1-kappa)/kappa))
eta_c_dyn = eta_c*(1 / (1 + exp(-(n_tc-30000)/(5000))))
P_c = mdot_beta*c_p*T_amb*1/eta_c_dyn*(p_im**((kappa-1)/kappa)-1)
P_f = 155

f = vertcat(u1_r,
     u2_r,
     1e-5 * 1/(1e5*xD1_scaled*J)*(P_t-P_c-P_f),
     1e-5 * R_l*T_ut/(150*V_ident)*(m_dot_exh-m_dot_egr_hp-mdot_turb)+0.001*u2)

g = vertcat(xA1, xA2)

sim = integrator('sim', 'idas', {'x': states, 'p': vertcat(controls, d1), 'z': z, 'ode': f, 'alg': g}, {'tf': Ts})

n_sim = 100
x_current = [1, 1, 0.1, 0.1]
for i in range(n_sim):
    print(x_current)
    x_current = sim(x0=x_current, p=[0, 0, 0, 2000])['xf']

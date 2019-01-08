
from casadi import *
from casadi.tools import *

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

def smoothen(x, eps=1e-4):
    return (sqrt(x**2 + eps) + x)/2

def smooth_fun(x, p, a):
    return a[0] + a[1] / (1 + exp(-(x-p[0])/p[1]))

def output(x, u):
    return vertcat(1000 * x['xD', 0] * x['xD', 1], x, u)

def output_N(x):
    return vertcat(1000 * x['xD', 0] * x['xD', 1], x)

def plot():
    plt.subplot(3, 1, 1)
    plt.plot(X[:, :2])
    plt.ylabel('controls')
    plt.subplot(3, 1, 2)
    plt.plot(X[:, 2:])
    plt.ylabel('states')
    plt.subplot(3, 1, 3)
    plt.plot(X[:, 2] * X[:, 3] * 1000)
    plt.plot(reference)
    plt.ylabel('output')

x = struct_symMX([
    entry('u', shape=2), entry('xD', shape=2)
])
u1 = x['u'][0]
u2 = x['u'][1]
xD1 = x['xD'][0]
xD2 = x['xD'][1]

x_ref = [50, 50, 1.14275, 1.53787]

d_diff_states = struct_symMX([
    entry('u', shape=2), entry('xD', shape=2)
])

u = struct_symMX([
    entry('u_r', shape=2)
])
u1_r = u['u_r'][0]
u2_r = u['u_r'][1]

alg_states = struct_symMX([
    entry('z', shape=2)
])
xA1 = alg_states['z'][0]
xA2 = alg_states['z'][1]

d1, d2 = (2000, 0)
y_ref = MX.sym('y_ref')

nx = x.size
nu = u.size
nz = alg_states.size

N = 20
Ts = 0.05

c1 = 25.3
c2 = 0.0034
c3 = 7.7e3
c4 = 0.6
c5 = 43.6
c6 = 9.2e-3
c7 = 3.6e3
c8 = 0.9

h_data = sio.loadmat('h_data.mat')
g_data = sio.loadmat('g_data.mat')
reference = np.ndarray.flatten(sio.loadmat('reference.mat')['reference'])

a = np.array([0.0, 1.0])
p = np.ndarray.flatten(h_data['p'])
g = np.ndarray.flatten(g_data['g'])
b = np.array([1.0, -1.0])

xA1_s = smoothen(xA1)
xA2_s = smoothen(xA2)

u1_star = smooth_fun(xD1*xD2+d2, p, a)*smooth_fun(u1, g, b)
u2_star = 1-(u2/100)

ode = vertcat(u1_r,
              u2_r,
              c1*(xA1_s**(1.5) - xA1_s**(1.25))*sqrt(smoothen(xA1_s**(-1.5) - xA1_s**(-1.75))) - c2*d1*xD2*(smoothen(xD1)**(1.29) - xD1),
              c5*xA1_s*(xA2_s**(1.5) - xA2_s**(1.25))*sqrt(smoothen(xA2_s**(-1.5) - xA2_s**(-1.75))) - c6*d1*xD1*(smoothen(xD2)**(1.29) - xD2))
            
alg = vertcat(-xD1*xD2 + c3/d1*sqrt(smoothen(xA1_s**(0.5) - xA1_s**(0.25)))*(xA1_s**(0.5) + c4*u1_star),
              -xD1*xD2 + c7/d1*xA1_s*sqrt(smoothen(xA2_s**(0.5) - xA2_s**(0.25)))*(xA2_s**(0.5) + c8*u2_star))

impl_dae = vertcat(d_diff_states - ode, alg)
jac_x = jacobian(impl_dae, x)
jac_d_x = jacobian(impl_dae, d_diff_states)
jac_u = jacobian(impl_dae, u)
jac_z = jacobian(impl_dae, alg_states)

inputs = [x, d_diff_states, u, alg_states]
engine_impl_dae_fun = Function('engine_impl_dae_fun', inputs, [impl_dae])
engine_impl_dae_fun_jac_x_xdot_z = Function('engine_impl_dae_fun_jac_x_xdot_z', inputs, [impl_dae, jac_x, jac_d_x, jac_z])
engine_impl_dae_jac_x_xdot_u_z = Function('engine_impl_dae_jac_x_xdot_u_z', inputs, [jac_x, jac_d_x, jac_u, jac_z])
# only needed for lifted IRK
engine_impl_dae_fun_jac_x_xdot_u_z = Function('engine_impl_dae_fun_jac_x_xdot_u_z', inputs, [impl_dae, jac_x, jac_d_x, jac_u, jac_z])
# objective residual
engine_ls_cost = Function('engine_ls_cost', [x,u], [output(x, u), jacobian(output(x, u), vertcat(u, x)).T])
engine_ls_cost_N = Function('engine_ls_cost_N', [x], [output_N(x), jacobian(output_N(x), x).T])

codegen_opts = {'mex': False, 'casadi_int': 'int', 'with_header': True}
for fun in [engine_impl_dae_fun, engine_impl_dae_fun_jac_x_xdot_z, engine_impl_dae_jac_x_xdot_u_z,
            engine_impl_dae_fun_jac_x_xdot_u_z, engine_ls_cost, engine_ls_cost_N]:
    fun.generate(fun.name(), codegen_opts)

sim = integrator('sim', 'collocation', {'x': x, 'p': u, 'z': alg_states, 'ode': ode, 'alg': alg},
                 {'tf': Ts, 'rootfinder': 'newton',  'number_of_finite_elements': 1, 'interpolation_order': 2})

V = struct_symMX([(
    entry('x', struct=x, repeat=N+1),
    entry('u', struct=u, repeat=N)
)])

constraints = struct_symMX([(
    entry('dynamics', struct=x, repeat=N)
)])

G = struct_MX(constraints)

x_current = [50, 50, 1.3244, 0.9568]
z_current = [1, 1]

# steady state
# x_ss, z_ss = [50, 50, 1.14275, 1.53787], [1.28976, 1.78264]

Q = np.eye(nx)
R = 1e-1*np.eye(nu)

objective = 0.0

for i in range(N):
    sim_out = sim(x0=V['x', i], z0=[1.28976, 1.78264], p=V['u', i])

    G['dynamics', i] = V['x', i+1] - sim_out['xf']

    objective += 100*(1000 * V['x', i, 'xD', 0] * V['x', i, 'xD', 1] - y_ref)**2
    objective += mtimes((V['x', i]-x_ref).T, mtimes(Q, (V['x', i])-x_ref))
    objective += mtimes(V['u', i].T, mtimes(R, V['u', i]))

objective += 100*(1000 * V['x', N, 'xD', 0] * V['x', N, 'xD', 1] - y_ref)**2
objective += mtimes((V['x', N]-x_ref).T, mtimes(Q, (V['x', N])-x_ref))

solver = nlpsol('solver', 'ipopt', {'x': V, 'f': objective, 'g': G, 'p': y_ref})

# bounds
lb, ub = V(-inf), V(+inf)

lb['x', :, 'u'] = 0.0
ub['x', :, 'u'] = 100.0
lb['x', :, 'xD'] = 0.5
ub['x', :, 'xD'] = repeated([1.757, 2.125])
lb['u', :, 'u_r'] = -10000.0
ub['u', :, 'u_r'] = +10000.0

# initial value
x_init = [50, 50, 1.3244, 0.9568]

# initial guess
x0 = V(0)

x0['x'] = repeated([50, 50, 1.3244, 0.9568])
x0['u'] = 0

X = [DM(x_init)]
U = []

for i in range(reference.size):

    lb['x', 0] = X[-1]
    ub['x', 0] = X[-1]

    solver_out = solver(x0=x0, lbx=lb, ubx=ub, lbg=0, ubg=0, p=reference[i])
    primal_solution = V(solver_out['x'])

    x_opt = primal_solution['x']
    u_opt = primal_solution['u']

    X.append(x_opt[1])
    U.append(u_opt[0])

    # shifting
    for i in range(N-1):
        x0['x', i] = primal_solution['x', i+1]
        x0['u', i] = primal_solution['u', i+1]
    x0['x', N-1] = primal_solution['x', N]
    x0['x', N] = primal_solution['x', N]

X = np.array(horzcat(*X).T)
U = np.array(horzcat(*U).T)

plt.ion()
plot()

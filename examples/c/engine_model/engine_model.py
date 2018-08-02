
from casadi import *
from casadi.tools import *

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

def smoothen(x, eps=1e-4):
    return (sqrt(x**2 + eps) + x)/2

def smooth_fun(x, p, a):
    return a[0] + a[1] / (1 + exp(-(x-p[0])/p[1]))

def plot():
    plt.subplot(3, 1, 1)
    plt.plot(X[:, :2])
    plt.subplot(3, 1, 2)
    plt.plot(X[:, 2:])
    plt.subplot(3, 1, 3)
    plt.plot(X[:, 2] * X[:, 3] * 1000)

class DebugCallback(Callback):
    def __init__(self, name, nx, ng, np=0, opts={}):

        self.nx = nx
        self.ng = ng
        self.np = np

        Callback.__init__(self)

        # Initialize internal objects
        self.construct(name, opts)

    def get_n_in(self): return nlpsol_n_out()
    def get_n_out(self): return 1
    def get_name_in(self, i): return nlpsol_out(i)
    def get_name_out(self, i): return "ret"

    def get_sparsity_in(self, i):
        n = nlpsol_out(i)
        if n=='f':
            return Sparsity.scalar()
        elif n in ('x', 'lam_x'):
            return Sparsity.dense(self.nx)
        elif n in ('g', 'lam_g'):
            return Sparsity.dense(self.ng)
        else:
            return Sparsity(0,0)
    
    def eval(self, arg):
        pass
        # Create dictionary

        # import ipdb; ipdb.set_trace()

        # darg = {}
        # print(arg[0])

diff_states = struct_symMX([
    entry('u', shape=2), entry('xD', shape=2)
])
u1 = diff_states['u'][0]
u2 = diff_states['u'][1]
xD1 = diff_states['xD'][0]
xD2 = diff_states['xD'][1]

d_diff_states = struct_symMX([
    entry('u', shape=2), entry('xD', shape=2)
])

controls = struct_symMX([
    entry('u_r', shape=2)
])
u1_r = controls['u_r'][0]
u2_r = controls['u_r'][1]

alg_states = struct_symMX([
    entry('z', shape=2)
])
xA1 = alg_states['z'][0]
xA2 = alg_states['z'][1]

d1, d2 = (2000, 0)
y_ref = MX.sym('y_ref')

nx = diff_states.size
nu = controls.size
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
jac_x = jacobian(impl_dae, diff_states)
jac_d_x = jacobian(impl_dae, d_diff_states)
jac_u = jacobian(impl_dae, controls)
jac_z = jacobian(impl_dae, alg_states)

inputs = [diff_states, d_diff_states, controls, alg_states]
engine_impl_dae_fun = Function('engine_impl_dae_fun', inputs, [impl_dae])
engine_impl_dae_fun_jac_x_xdot_z = Function('engine_impl_dae_fun_jac_x_xdot_z', inputs, [impl_dae, jac_x, jac_d_x, jac_z])
engine_impl_dae_jac_x_xdot_u_z = Function('engine_impl_dae_jac_x_xdot_u_z', inputs, [jac_x, jac_d_x, jac_u, jac_z])
# only needed for lifted IRK
engine_impl_dae_fun_jac_x_xdot_u_z = Function('engine_impl_dae_fun_jac_x_xdot_u_z', inputs, [impl_dae, jac_x, jac_d_x, jac_u, jac_z])

codegen_opts = {'mex': False, 'casadi_int': 'int', 'with_header': True}
for fun in [engine_impl_dae_fun, engine_impl_dae_fun_jac_x_xdot_z, engine_impl_dae_jac_x_xdot_u_z, engine_impl_dae_fun_jac_x_xdot_u_z]:
    fun.generate(fun.name(), codegen_opts)

sim = integrator('sim', 'collocation', {'x': diff_states, 'p': controls, 'z': alg_states, 'ode': ode, 'alg': alg},
                 {'tf': Ts, 'rootfinder': 'newton', 'number_of_finite_elements': 1, 'interpolation_order': 2})

V = struct_symMX([(
    entry('diff_states', struct=diff_states, repeat=N+1),
    entry('controls', struct=controls, repeat=N)
)])

constraints = struct_symMX([(
    entry('dynamics', struct=diff_states, repeat=N)
)])

G = struct_MX(constraints)

x_current = [50, 50, 1.3244, 0.9568]
z_current = [1, 1]

# steady state
# x_current, z_current = [50, 50, 1.14275, 1.53787], [1.28976, 1.78264]

Q = np.eye(nx)
R = 1e-1*np.eye(nu)

objective = 0.0

for i in range(N):
    sim_out = sim(x0=V['diff_states', i], z0=[1, 1], p=V['controls', i])

    G['dynamics', i] = V['diff_states', i+1] - sim_out['xf']

    objective += 100*(1000 * V['diff_states', i, 'xD', 0] * V['diff_states', i, 'xD', 1] - y_ref)**2
    objective += mtimes(V['diff_states', i].T, mtimes(Q, V['diff_states', i]))
    objective += mtimes(V['controls', i].T, mtimes(R, V['controls', i]))

objective += 100*(1000 * V['diff_states', N, 'xD', 0] * V['diff_states', N, 'xD', 1] - y_ref)**2
objective += mtimes(V['diff_states', N].T, mtimes(Q, V['diff_states', N]))

callback = DebugCallback('engine', V.size, G.size)
solver = nlpsol('solver', 'ipopt', {'x': V, 'f': objective, 'g': G, 'p': y_ref})

# bounds
lb, ub = V(-inf), V(+inf)

lb['diff_states', :, 'u'] = 0.0
ub['diff_states', :, 'u'] = 100.0
lb['diff_states', :, 'xD'] = 0.5
ub['diff_states', :, 'xD'] = repeated([1.757, 2.125])
lb['controls', :, 'u_r'] = -10000.0
ub['controls', :, 'u_r'] = +10000.0

# initial value
x_init = [50, 50, 1.3244, 0.9568]

# initial guess
x0 = V(0)

x0['diff_states'] = repeated([50, 50, 1.3244, 0.9568])
x0['controls'] = 0

X = [DM(x_init)]
U = []

for i in range(reference.size):

    lb['diff_states', 0] = X[-1]
    ub['diff_states', 0] = X[-1]

    solver_out = solver(x0=x0, lbx=lb, ubx=ub, lbg=0, ubg=0, p=reference[i])
    x0 = solver_out['x']

    x_opt = V(solver_out['x'])['diff_states']
    u_opt = V(solver_out['x'])['controls']

    X.append(x_opt[1])
    U.append(u_opt[0])

X = np.array(horzcat(*X).T)
U = np.array(horzcat(*U).T)

plt.ion()
plot()

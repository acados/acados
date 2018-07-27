from sys import path
path.append("/Users/robin/casadi-py27-np1.9.1-v2.4.3")
from casadi import *
from casadi.tools import *
import matplotlib.pyplot as plt
from time import sleep
import numpy as np
import numpy.random as random
import pickle

def pendulum_position(p,theta,l):
    return vertcat([p-l*sin(theta),l*cos(theta)])

def KKTnorm(df,dg,lam,g):
    KKT_residual = vertcat([df + mul(dg.T,lam),
                            g["initial_state"],
                            vertcat(g["dynamics"]),
                            max(0.,g["final_pendulum_position"]),
                            min(0.,lam["final_pendulum_position"]),
                            g["final_pendulum_position"]*lam["final_pendulum_position"]])
    return float(norm_2(KKT_residual))

def pendulum_final_position_constraint(p,theta,l):
    pendulum_final_position = pendulum_position(p,theta,l)
    return quad_form(pendulum_final_position - vertcat([l,l])) - radius2

states = struct_symMX(["p","pDot","theta","thetaDot"])
[p, pDot, theta, thetaDot] = states[...]
controls = struct_symMX(["F"])
[F] = controls[...]

# Problem parameters
M = 1.
m = 0.1
l = 0.8
g = 9.81
n_x = states.shape[0]
radius2 = 0.04

denominator = M + m - m*cos(theta)*cos(theta)
f_x = struct_MX(states)
f_x["p"] = pDot
f_x["pDot"] = (-m*l*sin(theta)*thetaDot*thetaDot + m*g*cos(theta)*sin(theta)+F)/denominator
f_x["theta"] = thetaDot
f_x["thetaDot"] = (-m*l*cos(theta)*sin(theta)*thetaDot*thetaDot + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator)

ode = MXFunction("ode",[states,controls],[f_x])
intg = simpleRK(ode,1,4)

# Optimization variables
n_sim = 20
W = struct_symMX([  (entry("states",struct=states,repeat=n_sim+1),
                    entry("controls",struct=controls,repeat=n_sim))])
constraints = struct_symMX([entry("initial_state",struct=states),
                            entry("dynamics",struct=states,repeat=n_sim),
                            # entry("final_velocity",shape=2),
                            entry("final_pendulum_position")])

# Optimization parameters
dt = 0.05
x0 = [0.,0.,pi,0.]
r_F = 1e-4
q_x = 1e-10

g = struct_MX(constraints)
g["initial_state"] = x0 - W["states",0]
g["final_pendulum_position"] = pendulum_final_position_constraint(W["states",n_sim,"p"],W["states",n_sim,"theta"],l)
obj = 0.
for i in range(n_sim):
    x_next = intg(x0=W["states",i],p=W["controls",i],h=dt)["xf"]
    g["dynamics",i] = x_next - W["states",i+1]
    # Control penalization
    obj += 0.5*r_F*quad_form(W["controls", i]) + 0.5*q_x*quad_form(W["states", i])

obj += 0.5*q_x*quad_form(W["states", n_sim])

nlp  = MXFunction("nlp",nlpIn(x=W),nlpOut(f=obj,g=g))

# Solve NLP with ipopt
tol = 1e-13
solver_options = {}
solver_options["linear_solver"] = "ma86"
solver_options["expand"] = True
solver_options["tol"] = tol
solver = NlpSolver("solver","ipopt",nlp,solver_options)
initial_guess = W(0.)
initial_guess["states",:] = repeated(x0)
old_sol = pickle.load(open("SCQPsol.p","r"))
initial_guess = old_sol[0]
lbg = constraints(0.)
lbg["final_pendulum_position"] = -inf
ubg = constraints(0.)

sol = solver({"x0":initial_guess,"lbg":lbg,"ubg":ubg})
prim_sol = sol["x"]
dual_sol = sol["lam_g"]

# Solve NLP with SQP
lam = struct_symMX(constraints)
L = MXFunction("L",[W,lam],[obj+mul(transpose(lam),g)])
exact_Hessian = MXFunction("exact_Hessian",[W,lam],[hessian(L([W,lam])[0],W)[0]])

g_fun = MXFunction("g_fun",[W],[g])
G = MXFunction("G",[W],[jacobian(g,W)])
Jc = jacobian(pendulum_position(W["states",n_sim,"p"],W["states",n_sim,"theta"],l),W)
Jc_function = MXFunction("Jc_function",[W],[Jc])
d2c = MXFunction("d2c",[W,lam],[lam["final_pendulum_position"]*2*mul(Jc.T,Jc)])
grad_obj = MXFunction("grad_obj",[W],[gradient(obj,W)])

# Hessian approximations
GN_Hessian = MXFunction("GN_Hessian",[W,lam],[hessian(obj,W)[0]])
SCQP_Hessian = MXFunction("SCQP_Hessian",[W,lam],[GN_Hessian([W,lam])[0] + d2c([W,lam])[0]])
hessians = [SCQP_Hessian]
max_iters = 100
j = 0
KKTtols = np.inf*np.ones([max_iters+1,3])
qpOASES_options = {"printLevel":"none","enableEqualities":True,"initialStatusBounds":"inactive","terminationTolerance":tol}
wk_close_to_optimum   = W(prim_sol +0.0001*random.randn(prim_sol.shape[0],1))
lamk_close_to_optimum = lam(dual_sol+0.0001*random.randn(dual_sol.shape[0],1))
for hessian_function in hessians:
    qpOASES_sparsity = {'h':hessian_function([W,lam])[0].sparsity(),'a':G([W])[0].sparsity()}
    qpsolver = QpSolver("qpsolver","qpoases",qpOASES_sparsity,qpOASES_options)
    wk = W(0)
    lamk = lam(0)
    k = 0
    KKTtol = KKTnorm(grad_obj([wk])[0],G([wk])[0],lamk,constraints(g_fun([wk])[0]))
    KKTtols[0,:] = KKTtol
    print "First KKTtol: ", KKTtol
    print "\t---" + hessian_function.getSanitizedName() + "---"
    print "i \t KKT norm \t ||w-w_star||"
    while KKTtol > tol and k < max_iters:
        [H] = hessian_function([wk,lamk])
        [A] = G([wk])
        [df] = grad_obj([wk])
        gg = constraints(g_fun([wk])[0])
        lba = constraints(-gg.cat)
        uba = constraints(-gg.cat)
        lba["final_pendulum_position"] = -inf
        dw = qpsolver(h=H,g=df,a=A,lba=lba,uba=uba)
        wk   = W(wk + dw["x"])
        lamk = lam(dw["lam_a"])
        KKTtol = KKTnorm(grad_obj([wk])[0],G([wk])[0],lamk,constraints(g_fun([wk])[0]))
        KKTtols[k+1,j] = KKTtol
        k += 1
        print '{0:3d}\t{1:3e}\t{2:6e}'.format(k, KKTtol, float(norm_2(wk-prim_sol)))
        # import ipdb; ipdb.set_trace()
    j += 1
    del(qpsolver)

ipopt_sol = W(prim_sol)
sqp_sol = W(wk)
sqp_lam = lam(lamk)
solutions = [ipopt_sol, sqp_sol]

plt.rc('font',**{'family':'serif','serif':'computer modern roman','size':18})
plt.rc('text',usetex=True)
# Plot KKT
plt.figure(3);plt.clf()
plt.semilogy(KKTtols[:,0],"k-")
plt.semilogy(KKTtols[:,1],"k--")
plt.semilogy(KKTtols[:,2],"k-.")
plt.legend(["EH","GGN","SCQP"],loc='best')
plt.title("$\| \mathrm{KKT~residual} \|$")
plt.ylim(1e-12,1e4)
plt.xlabel("iteration number")
plt.xlim(0,45)
plt.show()

plt.ion()
plt.figure(1);plt.clf()
plt.figure(2);plt.clf()
for k,sol in enumerate(solutions):
    plt.figure(1)
    plt.subplot(2,1,k+1,aspect=1.)
    for i in range(n_sim+1):
        p_i = sol["states",i,"p"]
        theta_i = sol["states",i,"theta"]
        pendulumXY = np.array([[p_i,0.],list(pendulum_position(p_i,theta_i,l))])
        color_i = 3*[0.9*(1-float(i)/n_sim)]
        plt.plot(pendulumXY[0,0],pendulumXY[0,1],"s",markersize=10.,color=color_i)
        plt.plot(pendulumXY[:,0],pendulumXY[:,1],"-",linewidth=2.,color=color_i)
        plt.plot(pendulumXY[1,0],pendulumXY[1,1],"o",markersize=5.,color=color_i)
    tt = np.linspace(0.,2*pi,50)
    plt.plot(sqrt(radius2)*cos(tt)+0.8,sqrt(radius2)*sin(tt)+l,'r-')
    plt.xlabel("X~[\mathrm{m}]")
    plt.ylabel("Y~[\mathrm{m}]")
    plt.ylim(-0.8,1)

    ax = plt.figure(2).gca()
    ax.set_aspect(0.008)
    plt.grid("on")
    plt.step(np.array(linspace(0,(n_sim-1)*dt,n_sim)),sol["controls"],'k-',where="post")
    plt.plot([(n_sim-1)*dt,n_sim*dt],[sol["controls",n_sim-1],sol["controls",n_sim-1]],'k-')
    plt.xlabel("$\mathrm{time}~[\mathrm{s}]$")
    plt.ylabel("$F~[\mathrm{N}]$")

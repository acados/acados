#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

# reference : "Towards Time-optimal Tunnel-following for Quadrotors", Jon Arrizabalaga et al.

import numpy as np
from matplotlib import pyplot, animation

from common import *
from sys_dynamics import SysDyn
from acados_template import latexify_plot
latexify_plot()

def animOptVars(misc_steps, traj_ST, traj_U):
    '''Plot data as animate matplotlib graph'''
    sysModel = SysDyn()
    zetaMx, _, _, _, _ = sysModel.SetupOde()

    # (nx , N+1 , mpc_iter)
    dim_st = np.shape(traj_ST)
    anim_running = True

    # Plot only the original track without repetition
    [_, xref_track, yref_track, zref_track] = getTrack()
    # subsample the solutions for faster animation

    # nx X N X (mpc_iter/freq)
    traj_ST = traj_ST[:, :, ::f_plot]
    # nu X N X (mpc_iter/freq)
    traj_U = traj_U[:, ::f_plot]
    # Contains miscelleanous info for plots (times, cost) 2x (mpc_iter/freq)
    misc_steps = misc_steps[:, ::f_plot]

    zetaC_hat = np.zeros((3, dim_st[1], dim_st[2]))
    zetaC_hat = zetaC_hat[:, :, ::f_plot]

    zetaEul = np.zeros((3, dim_st[1], dim_st[2]))
    zetaCEul = zetaEul[:, :, ::f_plot]

    list_kap = []
    list_tau = []
    kap, tau, et, en, eb, kapBar, tauBar, _, _, _ = projFrenSerretBasis(zetaMx[0])
    kap_fun, tau_fun, _, _, _ = evalFrenSerretBasis(zetaMx[0], kap, tau, et, en, eb)
    # kapBar_fun, tauBar_fun = evalFrenSerretDerv(zetaMx[0], kapBar, tauBar)

    for k in range(N+1):
        s_i = traj_ST[0, k, :]
        n_i = traj_ST[1, k, :]
        b_i = traj_ST[2, k, :]

        sDot_i = traj_ST[7, k, :]
        nDot_i = traj_ST[8, k, :]
        bDot_i = traj_ST[9, k, :]

        x_i, y_i, z_i = sysModel.Fren2CartT(zetaMx, s_i, n_i, b_i)

        zetaC_hat[0, k, : ] = np.ravel(DM2Arr(x_i))
        zetaC_hat[1, k, : ] = np.ravel(DM2Arr(y_i))
        zetaC_hat[2, k, : ] = np.ravel(DM2Arr(z_i))

        q1_i = traj_ST[3, k, :]
        q2_i = traj_ST[4, k, :]
        q3_i = traj_ST[5, k, :]
        q4_i = traj_ST[6, k, :]

    for i in range((s_ref.shape[0])):
       #et_fun1 = et_fun(s_ref[i]/2)
       #print("\t", gam4(s_ref[i]))
       list_kap.append(float(kap_fun(s_ref[i])))
       list_tau.append(float(tau_fun(s_ref[i])))

    def init():
        time_text.set_text('')

        return path, horizon, drone

    def onClick(event):
        nonlocal anim_running
        anim_running ^= True
        if anim_running:
            anim.event_source.stop()
            anim_running = False
        else:
            anim.event_source.start()
            anim_running = True

    def animate(iter):
        '''Update animation'''

        # update simulation time
        time_text.set_text(f'time = {misc_steps[0, iter]:.2f} s' )

        drone[0]._offsets3d = ([float(zetaC_hat[0, 0, iter])],
                       [float(zetaC_hat[1, 0, iter])],
                       [float(zetaC_hat[2, 0, iter])])

        horizon.set_data(zetaC_hat[0: 2, 1:, iter])
        horizon.set_3d_properties(zetaC_hat[2, 1:, iter])

        # Update state plots
        sAx.set_data(traj_ST[0,  0, : iter], misc_steps[0, :iter +1] )
        nAx.set_data(traj_ST[1,  0, : iter], misc_steps[0, :iter +1] )
        bAx.set_data(traj_ST[2,  0, : iter], misc_steps[0, :iter +1] )

        sDotAx.set_data(traj_ST[7,  0, : iter], misc_steps[0, :iter +1] )
        nDotAx.set_data(traj_ST[8,  0, : iter], misc_steps[0, :iter +1] )
        bDotAx.set_data(traj_ST[9,  0, : iter], misc_steps[0, :iter +1] )

        # Update control plot
        ohm1.set_data(traj_ST[-4, 0, :iter], misc_steps[0, :iter +1] )
        ohm2.set_data(traj_ST[-3, 0, :iter], misc_steps[0, :iter +1] )
        ohm3.set_data(traj_ST[-2, 0, :iter], misc_steps[0, :iter +1] )
        ohm4.set_data(traj_ST[-1, 0, :iter], misc_steps[0, :iter +1] )

        # # update cost plot
        # costAx.set_data(misc_steps[1, :iter], misc_steps[0, :iter +1] )

        return path, horizon, drone

    # Create a figure which occupies the full screen
#    et = np.append(np.asarray(et_fun(s_i[i])).flatten(), 0 )
#    en = np.append(np.asarray(en_fun(s_i[i])).flatten(), 0 )
#    eb = np.asarray([0,  0, 1])

    #Rf = np.vstack((et, en, eb))
    fig = pyplot.figure(figsize=(15, 10))

    # plot state on right (merge top and bottom right. i.e subplots 2, 3, 5, 6)
    ax3d = fig.add_subplot(3, 3, (5, 9), projection='3d')
    ax3d.azim = -25
    ax3d.elev = 15
    fig.add_axes(ax3d)

    # time field
    time_text = ax3d.text2D(0.02, 0.95, '', transform=ax3d.transAxes)

    # reference trajectory
    pyplot.plot(xref_track, yref_track, zref_track, linestyle='dashed', marker = 'x', c='cornflowerblue', dashes=(5, 15), alpha=0.3)
    # path
    path = ax3d.plot([], [], 'b', alpha=0.5, linewidth=0.5)[0]
    # horizon
    horizon, = ax3d.plot([], [],'x-g', alpha=0.5)

    cage_x = [-2.5, 2.5]
    cage_y = [-2.5, 2.5]
    cage_z = [0, 2]

    ax3d.set_aspect('auto')
    ax3d.set_xlim3d(left = cage_x[0], right = cage_x[1])
    ax3d.set_ylim3d(bottom = cage_y[0], top = cage_y[1])
    ax3d.set_zlim3d(bottom = cage_z[0], top = cage_z[1])

    # Single covering sphere around drone
    drone = [ None, None, None ]
    drone[0] = ax3d.scatter(x_i[0], y_i[0], z_i[0], s = np.pi * rob_rad**2 * sphere_scale, c='lightcoral', alpha=0.45)

    ax3d.set_xlabel('X (m)')
    ax3d.set_ylabel('Y (m)')
    ax3d.set_zlabel('Z (m)')

    # State zeta_frenet at (1,1) top left
    zetaF = fig.add_subplot(3, 3, 1)
    zetaF.set_ylim( np.amin(np.ravel(traj_ST[0 : 3, 0, :-2])) - 0.2,
                    np.amax(np.ravel(traj_ST[0 : 3, 0, :-2])) + 0.2)
    zetaF.set_xlim(0, np.amax(misc_steps[0, :]) + 0.2)
    zetaF.set_xlabel('time (s)')
    zetaF.set_ylabel("$p^{f}$")
    zetaF.set_yscale('symlog')
    sAx = zetaF.stairs([], [0], baseline=None,label="$s$ ($m$)", color="teal" )
    nAx = zetaF.stairs([], [0], baseline=None,label="$n$ ($m$)", color="lightcoral")
    bAx = zetaF.stairs([], [0], baseline=None,label="$b$ ($m$)", color="plum")
    zetaF.grid()
    zetaF.legend(loc='upper right')

    zetaDotF = fig.add_subplot(3, 3, 4)
    zetaDotF.set_ylim( np.amin(np.ravel(traj_ST[7 : 10, 0, :-2])) - 0.2,
                    np.amax(np.ravel(traj_ST[7 : 10, 0, :-2])) + 0.2)
    zetaDotF.set_xlim(0, np.amax(misc_steps[0, :]) + 0.2)
    zetaDotF.set_xlabel('time (s)')
    zetaDotF.set_ylabel("$\\dot{p}^{f}$")
    sDotAx = zetaDotF.stairs([], [0], baseline=None,label="$\\dot{s}$ ($m s^{-1}$)", color="teal" )
    nDotAx = zetaDotF.stairs([], [0], baseline=None,label="$\\dot{n}$ ($m s^{-1}$)", color="lightcoral")
    bDotAx = zetaDotF.stairs([], [0], baseline=None,label="$\\dot{b}$ ($m s^{-1}$)", color="plum")
    zetaDotF.grid()
    zetaDotF.legend(loc='upper right')


    # plot control u
    u = fig.add_subplot(3, 3, 2)
    ohm1 = u.stairs([], [0], baseline=None,label="$\\Omega_{1}$ ($rad s^{-1}$)", color="lightcoral" )
    ohm2 = u.stairs([], [0], baseline=None,label="$\\Omega_{2}$ ($rad s^{-1}$)", color="plum")
    ohm3 = u.stairs([], [0], baseline=None,label="$\\Omega_{3}$ ($rad s^{-1}$)", color="darkseagreen" )
    ohm4 = u.stairs([], [0], baseline=None,label="$\\Omega_{4}$ ($rad s^{-1}$)", color="lightsteelblue")

    u.set_ylim(np.amin(np.ravel(traj_ST[-4:, 0, :-2])) - 0.02,
               np.amax(np.ravel(traj_ST[-4:, 0, :-2])) + 0.02)
    u.set_xlim(0, np.amax(misc_steps[0, :]) + 0.2)
    u.set_xlabel('time (s)')
    u.set_ylabel('$\\zeta^{u}$')
    u.grid()
    u.legend(ncol=2, loc='upper right')

    # plot curve parametrization
    param = fig.add_subplot(3, 3, 3)
    tauAx = param.stairs(list_tau[:-1], s_ref, baseline=None,label="$\\tau$", color="coral" )
    kapAx = param.stairs(list_kap[:-1], s_ref, baseline=None,label="$\\kappa$", color="teal")
    param.set_xlim(0, np.amax(s_ref) + 0.2)
    param.set_ylim(   ymin = np.amin(np.ravel(list_tau[:] + list_kap[:]))*1.10 ,
                    ymax = np.amax(np.ravel(list_tau[:] + list_kap[:]))*1.10)
    param.set_xlabel('$s\,$(m)')
    param.legend(loc='upper right')
    param.grid()


    fig.canvas.mpl_connect('button_press_event', onClick)
    anim = animation.FuncAnimation(fig=fig, func=animate,
                                init_func=init,
                                frames=len(misc_steps[0, :]),
                                interval=refresh_ms,
                                repeat=True,
                                blit=False)

    fig.tight_layout()

    # avoid plotting when running on Travis
    if os.environ.get("ACADOS_ON_CI") is None:
        pyplot.show()
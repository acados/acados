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

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


def get_latex_plot_params():
    params = {'backend': 'ps',
            'text.latex.preamble': r"\usepackage{gensymb} \usepackage{amsmath}",
            'axes.labelsize': 12,
            'axes.titlesize': 12,
            'legend.fontsize': 12,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'text.usetex': True,
            'font.family': 'serif'
    }

    return params


def plot_chain_position_traj(simX, yPosWall=None):
    plt.figure()
    nx = simX.shape[1]
    N = simX.shape[0]
    M = int((nx/3 -1)/2)
    # plt.title('Chain position trajectory')

    for i in range(M+1):
        plt.subplot(M+1, 3, 3*i+1)
        plt.ylabel('x')
        plt.plot(simX[:, 3*i])
        plt.grid(True)

        plt.subplot(M+1, 3, 3*i+2)
        plt.ylabel('y')
        plt.plot(simX[:, 3*i+1])
        if not yPosWall == None:
            plt.plot(yPosWall*np.ones((N,)))
        plt.grid(True)

        plt.subplot(M+1, 3, 3*i+3)
        plt.ylabel('z')
        plt.plot(simX[:, 3*i+2])
        plt.grid(True)



def plot_chain_velocity_traj(simX):
    plt.figure()
    nx = simX.shape[1]
    M = int((nx/3 -1)/2)

    simX = simX[:, (M+1)*3:]

    for i in range(M):
        plt.subplot(M, 3, 3*i+1)
        plt.plot(simX[:, 3*i])
        plt.ylabel('vx')
        plt.grid(True)

        plt.subplot(M, 3, 3*i+2)
        plt.plot(simX[:, 3*i+1])
        plt.ylabel('vy')
        plt.grid(True)

        plt.subplot(M, 3, 3*i+3)
        plt.plot(simX[:, 3*i+2])
        plt.ylabel('vz')
        plt.grid(True)


def plot_chain_control_traj(simU):
    plt.figure()
    # plt.title('Chain control trajectory, velocities of last mass')
    simU = np.vstack((simU[0,:], simU))

    t = np.array(range(simU.shape[0]))
    plt.subplot(3, 1, 1)
    plt.step(t, simU[:,0])
    plt.ylabel('vx')
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.step(t, simU[:,1])
    plt.ylabel('vy')
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.step(t, simU[:,2])
    plt.ylabel('vz')
    plt.grid(True)



def plot_chain_position(x, xPosFirstMass):

    if len(x.shape) > 1:
        x = x.flatten()
    if len(xPosFirstMass.shape) > 1:
        xPosFirstMass = xPosFirstMass.flatten()

    nx = x.shape[0]
    M = int((nx/3 -1)/2)

    pos = x[:3*(M+1)]
    pos = np.hstack((xPosFirstMass, pos))  # append fixed mass
    pos_x = pos[::3]
    pos_y = pos[1::3]
    pos_z = pos[2::3]

    fig = plt.figure()
    plt.subplot(3,1,1)
    plt.plot(pos_x)
    plt.title('x position')
    plt.xlabel('mass index ')
    plt.ylabel('mass position ')
    plt.grid(True)
    
    plt.subplot(3,1,2)
    plt.plot(pos_y)
    plt.title('y position')
    plt.xlabel('mass index ')
    plt.ylabel('mass position ')
    plt.grid(True)
    
    plt.subplot(3,1,3)
    plt.plot(pos_z)
    plt.title('z position')
    plt.xlabel('mass index ')
    plt.ylabel('mass position ')
    plt.grid(True)


def plot_chain_position_3D(X, xPosFirstMass, XNames=None):
    """
    X can be either chain state, or tuple of chain states
    Xnames is a list of strings
    """

    if not isinstance(X, tuple):
        X = (X,)

    if XNames is None:
        XNames = []
        for i in range(len(X)):
            XNames += ['pos' + str(i + 1)]

    if not isinstance(XNames, list):
        XNames = [XNames]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(xPosFirstMass[0], xPosFirstMass[1], xPosFirstMass[2], 'rx')
    for i, x in enumerate(X):

        if len(x.shape) > 1:
            x = x.flatten()
        if len(xPosFirstMass.shape) > 1:
            xPosFirstMass = xPosFirstMass.flatten()
        
        nx = x.shape[0]
        M = int((nx/3 -1)/2)
        pos = x[:3*(M+1)]
        pos = np.hstack((xPosFirstMass, pos))  # append fixed mass
        pos_x = pos[::3]
        pos_y = pos[1::3]
        pos_z = pos[2::3]
        
        ax.plot(pos_x, pos_y, pos_z, '.-', label=XNames[i])
        
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.legend()


def get_plot_lims(a):
    
    a_min = np.amin(a)
    a_max = np.amax(a)
    # make sure limits are not equal to each other
    eps = 1e-12
    if np.abs(a_min - a_max) < eps:
        a_min -= 1e-3 
        a_max += 1e-3 

    return (a_min, a_max)
    

def animate_chain_position(simX, xPosFirstMass, Ts=0.1, yPosWall=None):
    '''
    Creates animation of the chain, where simX contains the state trajectory.
    dt defines the time gap (in seconds) between two succesive entries.
    '''
    
    # chain positions
    Nsim = simX.shape[0]
    nx = simX.shape[1]
    M = int((nx/3 -1)/2)
    pos = simX[:,:3*(M+1)]

    pos_x = np.hstack((xPosFirstMass[0] * np.ones((Nsim,1)), pos[:, ::3]))
    pos_y = np.hstack((xPosFirstMass[1] * np.ones((Nsim,1)), pos[:, 1::3]))
    pos_z = np.hstack((xPosFirstMass[2] * np.ones((Nsim,1)), pos[:, 2::3]))


    # limits in all three dimensions
    
    # ylim_x = (np.amin( pos_x), np.amax( pos_x))
    # ylim_y = (np.amin( pos_y), np.amax( pos_y))
    # ylim_z = (np.amin( pos_z), np.amax( pos_z))
    # eps = 1e-12
    # if np.abs(ylim_x[0] - ylim_x[1]) < eps:
    #     ylim_x[0] += 1e-3 
    #     ylim_x[0] += 1e-3 

    ylim_x = get_plot_lims(pos_x)
    ylim_y = get_plot_lims(pos_y)
    if yPosWall is not None:
        ylim_y = (min(ylim_y[0], yPosWall) - 0.1, ylim_y[1])
    ylim_z = get_plot_lims(pos_z)

    fig = plt.figure()    
    ax1 = fig.add_subplot(311, autoscale_on=False, xlim=(0,M+2), ylim=ylim_x)
    plt.grid(True)
    ax2 = fig.add_subplot(312, autoscale_on=False, xlim=(0,M+2), ylim=ylim_y)
    plt.grid(True)
    ax3 = fig.add_subplot(313, autoscale_on=False, xlim=(0,M+2), ylim=ylim_z)
    plt.grid(True)

    ax1.set_ylabel('x')
    ax2.set_ylabel('y')
    ax3.set_ylabel('z')

    # ax.set_aspect('equal')
    # ax.axis('off')

    # create empty plot
    line1, = ax1.plot([], [], '.-')
    line2, = ax2.plot([], [], '.-')
    line3, = ax3.plot([], [], '.-')
    
    lines = [line1, line2, line3]
    
    if yPosWall is not None:
        ax2.plot(yPosWall*np.ones((Nsim,)))
    
        
    def init():
        # placeholder for data
        lines = [line1, line2, line3]
        for line in lines:
            line.set_data([],[])

        # lines[0].set_data(list(range(M+2)), pos_x[0,:])
        # lines[1].set_data(list(range(M+2)), pos_y[0,:])
        # lines[2].set_data(list(range(M+2)), pos_z[0,:])
        return lines


    def animate(i):

        lines[0].set_data(list(range(M+2)), pos_x[i,:])
        lines[1].set_data(list(range(M+2)), pos_y[i,:])
        lines[2].set_data(list(range(M+2)), pos_z[i,:])

        return lines

    ani = animation.FuncAnimation(fig, animate, Nsim,
                                  interval=Ts*1000, repeat_delay=500,
                                  blit=True, init_func=init)
    plt.show()
    return ani


def animate_chain_position_3D(simX, xPosFirstMass, Ts=0.1):
    '''
    Create 3D animation of the chain, where simX contains the state trajectory.
    dt defines the time gap (in seconds) between two succesive entries.
    '''
    
    # chain positions
    Nsim = simX.shape[0]
    nx = simX.shape[1]
    M = int((nx/3 -1)/2)
    pos = simX[:,:3*(M+1)]
    # import pdb; pdb.set_trace()
    pos_x = np.hstack((xPosFirstMass[0] * np.ones((Nsim,1)) , pos[:, ::3]))
    pos_y = np.hstack((xPosFirstMass[1] * np.ones((Nsim,1)),  pos[:, 1::3]))
    pos_z = np.hstack((xPosFirstMass[2] * np.ones((Nsim,1)),  pos[:, 2::3]))


    xlim = get_plot_lims(pos_x)
    ylim = get_plot_lims(pos_y)
    zlim = get_plot_lims(pos_z)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d', autoscale_on=False, xlim=xlim, ylim=ylim, zlim=zlim)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # ax.set_aspect('equal')
    # ax.axis('off')

    # create empty plot
    # line, = ax.plot([], [], [], '.-')
    line, = ax.plot(pos_x[0,:], pos_y[1,:], pos_z[2,:], '.-')

    def init():
        # placeholder for data
        # line.set_data([], [])
        # line.set_3d_properties([])
        return line,


    def animate(i):

        line.set_data(pos_x[i,:], pos_y[i,:])
        # line.set_data(pos_x[i,:], pos_y[i,:], pos_z[i,:])
        line.set_3d_properties(pos_z[i,:])
        return line,

    ani = animation.FuncAnimation(fig, animate, Nsim,
                                  interval=Ts*1000, repeat_delay=500,
                                  blit=True,)# init_func=init)
    plt.show()
    return ani
    

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
import casadi as ca
import os
from pathlib import Path
from typing import Union

'''Global variables'''

track="trefoil_track.txt"

def getTrack():
    track_file = os.path.join(str(Path(__file__).parent), "tracks/", track)
    array=np.loadtxt(track_file, skiprows=1)
    sref = array[1:,0]
    xref = array[1:,1]
    yref = array[1:,2]
    zref = array[1:,3]

    return sref, xref, yref, zref

[s_ref, x_ref, y_ref, z_ref] = getTrack()

length = len(s_ref)
pathlength = s_ref[-1]

# CrazyFlie 2.1 physical parameters
g0  = 9.80665       # [m.s^2] gravitational accerelation
mq  = 31e-3         # [kg] total mass (with Lighthouse deck)
Ix = 1.395e-5       # [kg.m^2] Inertial moment around x-axis
Iy = 1.395e-5       # [kg.m^2] Inertial moment around y-axis
Iz = 2.173e-5       # [kg.m^2] Inertia moment around z-axis
Cd  = 7.9379e-06    # [N/krpm^2] Drag coefficient
Ct  = 3.25e-4       # [N/krpm^2] Thrust coefficient
dq  = 92e-3         # [m] distance between motors' center
l   = dq/2          # [m] distance between motors' center and the axis of rotation

# timing parameters
T_del = 0.02               # time between steps in seconds
N = 50                     # number of shooting nodes
Tf = N * T_del

Tsim = 45
Nsim = int(Tsim * N / Tf)

U_MAX = 22                                              # [krpm]
U_HOV = int(np.sqrt(.25 * 1e6* mq * g0 /Ct)) /1000      #[krpm]
U_REF = np.array([U_HOV, U_HOV, U_HOV, U_HOV])

# State
n_states = 20

init_zeta = np.array([0.05, 0, 0,       # s,  n,  b
                      1, 0, 0, 0,       # qw, qx, qy, qz,
                      .02, 0, 0,        # sdot, ndot, bdot
                      0, 0, 0,          # ohmr,  ohmp,  ohmy
                      0, 0, 0,          # vx, vy, vz
                      U_HOV, U_HOV, U_HOV, U_HOV ])     # ohm1, ohm2, ohm3, ohm4

rob_rad = 0.04                           # radius of the drone covering sphere

# Control
n_controls = 4

# Weights & Tracking reference
S_REF = 0.1875
S_MAX = 5.9
                                          # State weights on
Q = np.diag([1, 1e-1, 1e-1,               # frenet position
             1e-5, 1e-5, 1e-5, 1e-5,      # quaternion
             1e-5, 1e-5, 1e-5,            # frenet velocity
             1e-5, 1e-5, 1e-5,            # drone angular velocity
             1e-5, 1e-5, 1e-5,            # cartesian velocity
             1e-8, 1e-8, 1e-8, 1e-8])     # rotor angular velocity

                                          # Terminal state weights on
Qn = np.diag([10, 1e-3, 1e-2,             # frenet position
             1e-5, 1e-5, 1e-5, 1e-5,      # quaternion
             1e-5, 1e-5, 1e-5,            # frenet velocity
             1e-5, 1e-5, 1e-5,            # drone angular velocity
             1e-5, 1e-5, 1e-5,            # cartesian velocity
             1e-8, 1e-8, 1e-8, 1e-8])     # rotor angular velocity

R = np.diag([1e-5, 1e-5, 1e-5, 1e-5])

#  MatPlotLib animation settings
Tstart_offset = 0
f_plot = 10
refresh_ms = 10
sphere_scale = 20000 #TODO make dependant on map size. (10000/ 20 obst)
z_const = 0.1

''' Helper functions'''

def DM2Arr(dm):
    return np.array(dm.full(), dtype=object)

def quat2rpy(qoid):
    ''' qoid -> [qw, qx, qy, qz]
        returns euler angles in degrees
        reference math3d.h crazyflie-firmware'''

    r	  =  ca.atan2( 2 * (qoid[0]*qoid[1] + qoid[2]*qoid[3]), 1 - 2 * (qoid[1]**2 + qoid[2]**2 ))
    p     =  ca.asin( 2 *  (qoid[0]*qoid[2] - qoid[1]*qoid[3]))
    y	  =  ca.atan2( 2 * (qoid[0]*qoid[3] + qoid[1]*qoid[2]), 1 - 2 * (qoid[2]**2 + qoid[3]**2 ))

    r_d = r * 180 / np.pi          # roll in degrees
    p_d = p * 180 / np.pi          # pitch in degrees
    y_d = y * 180 / np.pi          # yaw in degrees

    return [r_d, p_d, y_d]


def InterpolLuT(s: Union[ca.MX, float]):
    ''' 3rd bspline interpolation of curve x, y, zeta based on longitudinal progress (s)
    <-- xref, yref, zref : position reference curve interpol function'''

    x_ref_curve = ca.interpolant("x_ref", "bspline", [s_ref], x_ref)
    y_ref_curve = ca.interpolant("y_ref", "bspline", [s_ref], y_ref)
    z_ref_curve = ca.interpolant("z_ref", "bspline", [s_ref], z_ref)

    return x_ref_curve(s), y_ref_curve(s), z_ref_curve(s)

def projFrenSerretBasis(s: Union[ca.MX, float]):
    ''' project to the Frenet Serret space
    <-- kap_impl : Curvature
    <-- tau_impl : torsion
    <-- dGamma_ds, d2Gamma_ds2, d3Gamma_ds3 : First 3 derivates of the curve w.r.t arc length
    <-- et_MX, en_MX, eb_MX : tangent , normal , binormal unit vectors'''

    InterOpts = {'degree': [5]}
    y_ref_MX = ca.interpolant("y_ref", "bspline", [s_ref], y_ref, InterOpts)
    z_ref_MX = ca.interpolant("z_ref", "bspline", [s_ref], z_ref, InterOpts)
    x_ref_MX = ca.interpolant("x_ref", "bspline", [s_ref], x_ref, InterOpts)


    [d2GammaX_ds2, dGammaX_ds] = ca.hessian(x_ref_MX(s), s)
    [d2GammaY_ds2, dGammaY_ds] = ca.hessian(y_ref_MX(s), s)
    [d2GammaZ_ds2, dGammaZ_ds] = ca.hessian(z_ref_MX(s), s)

    [d4GammaX_ds4, d3GammaX_ds3] = ca.hessian(d2GammaX_ds2, s)
    [d4GammaY_ds4, d3GammaY_ds3] = ca.hessian(d2GammaY_ds2, s)
    [d4GammaZ_ds4, d3GammaZ_ds3] = ca.hessian(d2GammaZ_ds2, s)

    dGamma_ds = ca.vertcat(dGammaX_ds, dGammaY_ds, dGammaZ_ds)
    d2Gamma_ds2 = ca.vertcat(d2GammaX_ds2, d2GammaY_ds2, d2GammaZ_ds2)
    d3Gamma_ds3 = ca.vertcat(d3GammaX_ds3, d3GammaY_ds3, d3GammaZ_ds3)
    d4Gamma_ds4 = ca.vertcat(d4GammaX_ds4, d4GammaY_ds4, d4GammaZ_ds4)

    kap = ca.norm_2(d2Gamma_ds2)
    kapBar_MX = ca.jacobian(kap, s)

    et_MX = dGamma_ds
    en_MX = d2Gamma_ds2/ kap
    eb_MX = ca.cross(et_MX, en_MX)

    etBar_MX = d2Gamma_ds2
    enBar_MX = (1/ kap**3) * (kap**2 * d3Gamma_ds3 - d2Gamma_ds2 @ d3Gamma_ds3.T @ d2Gamma_ds2)
    ebBar_MX = (1/ kap) * ca.cross(dGamma_ds, d3Gamma_ds3)
    tau_MX =  ca.dot(en_MX, ebBar_MX)
    tauBar_MX =  ca.jacobian(tau_MX, s)

    return kap, tau_MX, et_MX, en_MX, eb_MX, kapBar_MX, tauBar_MX, etBar_MX, enBar_MX, ebBar_MX

def evalFrenSerretBasis(s: ca.MX, kap_MX, tau_MX, et_MX, en_MX, eb_MX):
    '''evaluation functions for curve in Frenet Serret space
    <-- kap_fun, tau_fun : Curvature, torsion CasADi functions
    <-- et_fun, en_fun, eb_fun : tangent , normal , binormal CasADi functions '''

    tau_fun = ca.Function('tau', [s], [tau_MX])
    kap_fun = ca.Function('kap', [s], [kap_MX])
    et_fun = ca.Function('et', [s], [et_MX])
    en_fun = ca.Function('en', [s], [en_MX])
    eb_fun = ca.Function('en', [s], [eb_MX])

    return kap_fun, tau_fun, et_fun, en_fun, eb_fun

def evalFrenSerretDerv(s: ca.MX, kapBar, tauBar):
    '''evaluation functions for curve in Frenet Serret space
    <-- kap_fun, tau_fun : Curvature, torsion CasADi functions
    <-- et_fun, en_fun, eb_fun : tangent , normal , binormal CasADi functions '''

    tauBar_fun = ca.Function('tau', [s], [kapBar])
    kapBar_fun = ca.Function('kap', [s], [tauBar])

    return kapBar_fun, tauBar_fun

